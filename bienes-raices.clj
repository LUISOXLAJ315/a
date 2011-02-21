(require '[clojure.string :as s]
         '[incanter.stats :as stats]
         '[clojure.contrib.generic.math-functions :as m]
         '[clojure.contrib.combinatorics :as c])

;;;; Conversion

(defn get-rows [s]
  (let [[sep & ss] (.split s "\n")]
    (take-nth 2 (partition-by #{sep} ss))))

;; El archivo pdf primero debe ser convertido a HTML usando pdftohtml
;; y luego hay que sacarle la parte HTML usando:
;; sed -ie 's/<br>//g' archivo.html
;; sed -ie 's/&nbsp;/ /g' archivo.html
;; vim archivo.html
(defn convert [input-filename output-filename]
  (binding [*out* (java.io.FileWriter. output-filename)]
    (println "Numero comunero;Nombres;Apellido paterno;Apellido materno")
    (doseq [[num-comunero nombres apellidos] (get-rows (slurp input-filename))
            :let [apellidos (when apellidos
                              (s/replace apellidos #"\ +" " "))]]
      (println (apply str (interpose ";" (concat [num-comunero nombres]
                                                 (when apellidos
                                                   (.split apellidos " ")))))))))


;;;; Probabilities

(defn simulate-prob-n [pred csv-name n & {:keys [eps n-success] :or {eps 1E-6 n-success 5}}]
  (let [lns (vec (for [line (rest (.split (slurp csv-name) "\n"))]
                   (first (drop 2 (.split line ";")))))
        gen-group #(stats/sample lns :size n :replacement false)
        gen-1000-gs #(repeatedly 1000 gen-group)]
    (loop [last-p 0
           i 0
           j 0
           gs (gen-1000-gs)]
      (let [p (/ (count (filter (partial apply pred) gs)) 1000.0)
            sum-p (/ (+ (* last-p n) p) (inc n))
            j (if (m/approx= sum-p last-p eps) (inc j) 0)]
        (if (>= j n-success)
          [p i]
          (recur sum-p (inc i) j (gen-1000-gs)))))))

(defn ningun-apellido-en-comun [& as]
  (not-any? (partial apply =) (partition 2 1 (sort as))))

(defn tree-prob-n [csv-name n]
  (let [lns (for [line (rest (.split (slurp csv-name) "\n"))]
              (first (drop 2 (.split line ";"))))]
    ))

(defn binomial-coefficient [n k]
  (if (> k n) 0 (/ (apply * (range (inc (- n k)) (inc n)))
                   (apply * (range 1 (inc k))))))

(defn comb-prob-n [csv-name n]
  (let [lns (for [line (rest (.split (slurp csv-name) "\n"))]
              (first (drop 2 (.split line ";"))))
        ;; for each lastname we can compute how many people have it
        freqs (frequencies lns)
        P (binomial-coefficient (count lns) n)
        sum (reduce + (for [[_ m] freqs
                            i (range 2 (inc n))
                            :let [in-group (binomial-coefficient m i)
                                  out-group (binomial-coefficient (- P m) (- n i))]]
                        (* in-group out-group)))]
    ;; FIXME this is wrong. it overestimates the probability because of overlapping between groups
    (/ sum P)))

(defn exhaustive-prob-n [pred csv-name n]
  (let [lns (for [line (rest (.split (slurp csv-name) "\n"))]
              (first (drop 2 (.split line ";"))))
        conteo (for [as (c/combinations lns n)]
                 (if (apply pred as) 1 0))]
    (/ (apply + conteo) (count conteo))))


(defn write-prob-table [input-filenames output-filename n]
  (binding [*out* (java.io.FileWriter. output-filename)]
    (println "Comunidad;Prob ningun apellido en comun")
    (doseq [i input-filenames
            :let [name (s/capitalize (first (.split (slurp i) "\n")))]]
      (println (str name ";" (float (prob-n i n)))))))

(defn file-list [start end error-filename]
  (let [erroneous (set (map #(Integer/parseInt %) (.split (slurp error-filename) "\n")))]
    (for [i (range start end) :when (not (erroneous i))]
      (str i ".txt"))))


;;;; Combinaciones no unicas tipo hermano y tipo dobles primos

(defn- print-p-and-cummulative-p [freqs & [get-surnames]]
  ;; freqs, que se obtuvo con la funcion frequencies, es una especie de diccionario en el que
  ;; se encuentran apareados los apellidos con la cantidad de veces que este fue encontrado
  ;; es decir freqs es algo como {["PATERNO1" "MATERNO1"] 5, ["PATERNO2" "MATERNO2"] 3, ...}
  ;; siendo ["PATERNO1" "MATERNO1"] una llave (key) y 5 el valor (val)
  (let [total (apply + (vals freqs)) ;; (vals freqs) nos devuelve los valores del diccionario,
        ;; que al sumarlos todos nos da la cantidad de personas en la comunidad
        cum (atom 0) ;; un atom es un valor que puede mutar (en Clojure la gran mayoria de los
        ;; valores son immutables) y que usaremos para ir guardando la frec rel acumulada
        p-and-cumulative-p (for [[lns freq] (sort-by val > freqs) ;; primero ordenamos las
                                 ;; frecuencias de mayor a menor
                                 ;; en lns estaran los apellidos y en freq la frecuencia
                                 :let [p (/ freq total)]] ;; calculamos el porcentaje que esa
                             ;; frecuencia representa de la poblacion total
                             ;; y luego devolvemos los apellidos sin modificar + la frecuencia +
                             ;; el porcentaje + el porcentaje acumulado
                             ;; swap! es la funcion que se usa para modificar el atom cum
                             [lns [freq p (swap! cum + p)]])
        ;; en otras palabras, en p-and-cumulative-p tenemos lo mismo que en freqs,
        ;; pero ordenado y con el porcentaje y porcentaje acumulado
        ;; en la sgte linea filtramos de freqs todos los apellidos que tengan frecuencia 1
        combinaciones-no-unicas (filter (fn [[_ freq]] (> freq 1)) freqs)
        ;; finalmente calculamos la suma de los pares de apellidos repetidos
        sum (apply + (map val combinaciones-no-unicas))]
    ;; lo que imprimimos a continuacion es lo que va a los *-h.csv y *-p.csv
    (println "Paterno;Materno;Frecuencia;Frecuencia relativa;Frecuencia relativa acumulada")
    (doseq [[lns [freq p cum]] p-and-cumulative-p ;; sencillamente sacamos la informacion de
            ;; p-and-cumulative-p
            ;; la sgte linea es un truco para sacar los apellidos en el caso de dobles primos
            :let [[ln1 ln2] ((or get-surnames identity) lns)]]
      ;; en la sgte linea imprimimos. p y cum los convertimos a numeros decimales,
      ;; pq estan como fracciones (que es una gran gracia de Clojure, pq los numeros decimales son
      ;; imprecisos)... podriamos imprimir la fraccion tb, pero es mas facil leer numeros decimales
      (println (str ln1 ";" ln2 ";" freq ";" (float p) ";" (float cum))))
    ;; lo sgte es el valor que devuelve la funcion, que nos serviria para hacer los calculos
    ;; de la tabla mas general que me mostrabas en el dibujo
    [sum (/ sum total)]))

(defn sibling-like [csv-name]
  ;; las dos primeras lineas leen el csv y sacan todos los apellidos,
  ;; filtrando las personas que no tienen dos apellidos
  ;; lnss es por ende una lista de pares de apellidos: [["PATERNO1" "MATERNO1"] ...]
  (let [lnss (filter next (for [line (rest (.split (slurp csv-name) "\n"))]
                            (filter seq (take 2 (drop 2 (.split line ";"))))))]
    ;; luego la funcion frequencies calcula cuantas veces se encuentra cada par de apellidos
    ;; en la lista. es importante que para eso usa la funcion =, de modo que los apellidos
    ;; tienen que estar en el mismo orden para que los cuente como el mismo
    (print-p-and-cummulative-p (frequencies lnss))))

(defn get-surnames-from-map [m]
  (mapcat (fn [[k v]] (repeat v k)) m))

(defn to-map [s]
  (into {} (for [x (distinct s)]
             [x (count (filter #{x} s))])))

(defn double-cousins-like [csv-name]
  ;; aqui es lo mismo de antes, pero hacemos una magia para que no importe el orden.
  ;; la magia es que convertimos los pares de apellidos ["PATERNO1" "MATERNO1"] a diccionarios
  ;; {"PATERNO1" 1, "MATERNO1" 1} y en los diccionarios no importa el orden cuando uno los
  ;; compara con =
  (let [lnss (map to-map (filter next (for [line (rest (.split (slurp csv-name) "\n"))]
                                        (filter seq (take 2 (drop 2 (.split line ";")))))))]
    (print-p-and-cummulative-p (frequencies lnss) get-surnames-from-map)))


(defn batch [start end error-filename & [input-prefix output-prefix]]
  (doseq [i (map #(apply str (drop-last 4 %)) (file-list start end error-filename))
          :let [in (str input-prefix i ".csv")]]
    (binding [*out* (java.io.FileWriter. (str output-prefix i "-h.csv"))]
      (sibling-like in))
    (binding [*out* (java.io.FileWriter. (str output-prefix i "-p.csv"))]
      (double-cousins-like in))))
