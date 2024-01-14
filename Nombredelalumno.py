from itertools import combinations
from collections import defaultdict
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import gudhi as gd 
from scipy.spatial.distance import pdist, squareform

#################################################
# Practicas 1 y 3
class ComplejoSimplicial:
    def __init__(self):
        """
        Inicializa un complejo simplicial vacío.
        """
        self.simplices = []  # Lista para almacenar los símplices y sus valores de filtración
        self.grafo = defaultdict(set)

    def añadir_simplex(self, simplex, filtracion):
        """
        Añade un simplex al complejo simplicial y lo ordena.
        Si el simplex es una arista, se añade al grafo para el cálculo de las componentes conexas.

        :param simplex: Un conjunto que representa un simplex (e.g., {1, 2, 3}).
        :param filtracion: Valor flotante asociado al simplex.
        """
        self.simplices.append((simplex, filtracion))
        self.simplices.sort(key=lambda x: (x[1], len(x[0])))
        # Agregar aristas al grafo para calcular componentes conexas
        if len(simplex) == 2:
            u, v = simplex
            self.grafo[u].add(v)
            self.grafo[v].add(u)

    def calcular_dimension(self):
        """
        Calcula la dimensión del complejo simplicial, basada en la dimensión de sus símplices.

        :return: La dimensión del complejo simplicial.
        """
        return max(len(simplex) - 1 for simplex, _ in self.simplices)
    

    def obtener_todas_las_caras(self):
        """
        Calcula el conjunto de todas las caras de cada símplice en el complejo simplicial.

        Un "símplice" en un complejo simplicial es una estructura que generaliza la noción de un triángulo o tetraedro a cualquier número de dimensiones. 
        Las "caras" de un símplice son todos los subconjuntos posibles de sus vértices que también forman símplices. 
        Por ejemplo, las caras de un triángulo (2-símplice) incluyen sus tres vértices (0-símplices), sus tres aristas (1-símplices), y el triángulo en sí.

        :return: Un conjunto que contiene todas las caras de todos los símplices en el complejo simplicial.
        """
        todas_las_caras = set()
        for simplex, _ in self.simplices:
            # Itera sobre cada símplice en el complejo simplicial.
            for i in range(len(simplex)):
                # Genera todas las combinaciones posibles de los vértices del símplice.
                # Estas combinaciones representan todas las caras posibles del símplice.
                for cara in combinations(simplex, i + 1):
                    # Agrega cada cara (como un conjunto inmutable 'frozenset') al conjunto de todas las caras.
                    # Usamos 'frozenset' porque las caras son conjuntos de vértices que no deben modificarse.
                    todas_las_caras.add(frozenset(cara))
        return todas_las_caras


    def caras_por_dimension(self, dimension):
        """
        Calcula y devuelve el conjunto de todas las caras de una dimensión específica en el complejo simplicial.

        En matemáticas, una "cara" de un símplice se refiere a cualquier símplice que es parte de él, incluyendo a sí mismo. 
        Por ejemplo, en un triángulo (2-símplice), las caras de dimensión 0 son sus vértices, las caras de dimensión 1 son sus aristas, 
        y una cara de dimensión 2 es el triángulo completo.

        :param dimension: La dimensión de las caras que se desea calcular. La dimensión se refiere al número de dimensiones 
        que tiene el símplice. Por ejemplo, 0 para puntos, 1 para líneas, 2 para triángulos, etc.

        :return: Un conjunto que contiene todas las caras del complejo simplicial que tienen exactamente la dimensión especificada.
        """
        return {cara for cara in self.obtener_todas_las_caras() if len(cara) - 1 == dimension}
        # Aquí se utiliza una comprensión de conjunto para filtrar y obtener solo las caras de la dimensión deseada.
        # 'self.obtener_todas_las_caras()' llama a la función que genera todas las caras del complejo.
        # 'len(cara) - 1' calcula la dimensión de cada cara. En matemáticas, la dimensión de un símplice se define como 
        # el número de vértices que tiene menos uno. Por ejemplo, un punto (0-símplice) tiene 1 vértice, 
        # una línea (1-símplice) tiene 2 vértices, un triángulo (2-símplice) tiene 3 vértices, y así sucesivamente.
        # Por lo tanto, si 'len(cara) - 1' es igual a la 'dimension' dada, esa cara es de la dimensión que estamos buscando 
        # y se incluye en el conjunto resultante.


    def calcular_estrella(self, simplex_buscado):
        """
        Calcula la "estrella" de un símplice en un complejo simplicial.

        En topología, la "estrella" de un símplice en un complejo simplicial se define como el conjunto de todos los símplices 
        que contienen al símplice dado como una cara. La estrella de un símplice nos da una idea de cómo está conectado este 
        símplice con otros símplices más grandes en el complejo.

        :param simplex_buscado: El símplice para el cual se calculará la estrella. Un símplice es un conjunto de vértices, 
        y aquí se espera que 'simplex_buscado' sea una colección de tales vértices.

        :return: Un conjunto de símplices que forman la estrella del símplice dado. Cada símplice en esta estrella contiene 
        al 'simplex_buscado' como una subestructura.

        """
        simplex_buscado = frozenset(simplex_buscado)
        # Se convierte el símplice buscado en un 'frozenset' para garantizar que sea un conjunto inmutable de vértices.

        estrella = set()
        for simplex, _ in self.simplices:
            # Itera sobre todos los símplices en el complejo simplicial.
            if simplex_buscado.issubset(simplex):
                # Verifica si el 'simplex_buscado' es una cara del símplice actual.
                # Si es así, el símplice actual es parte de la estrella del 'simplex_buscado'.
                estrella.add(frozenset(simplex))
                # Agrega el símplice a la estrella. Se usa 'frozenset' para mantener la inmutabilidad.

        return estrella
        # Devuelve el conjunto de símplices que forman la estrella del 'simplex_buscado'.

    
    def calcular_link(self, simplex_buscado):
        """
        Calcula el link de un símplice dado.

        El link de un símplice τ en un complejo simplicial K, Lk(τ), es el conjunto de todas
        las caras de los símplices en la estrella de τ que no tienen vértices en común con τ.

        :param simplex_buscado: El símplice para el cual se calculará el link.
        :return: Un conjunto con todos los símplices que forman el link del símplice dado.
        """
        simplex_buscado = frozenset(simplex_buscado)
        estrella = self.calcular_estrella(simplex_buscado)
        link = set()

        # Para cada símplice en la estrella, obtener todas sus caras que no son caras de τ
        for simplex in estrella:
            # Comprobar si el simplex actual es una cara de τ; si lo es, no pertenece al link
            if simplex_buscado.issubset(simplex) and simplex != simplex_buscado:
                for i in range(len(simplex)):
                    for cara in combinations(simplex, i + 1):
                        cara_set = frozenset(cara)
                        # Si la cara no contiene al simplex_buscado y no es el simplex_buscado mismo,
                        # entonces es parte del link
                        if not simplex_buscado.issubset(cara_set) and cara_set != simplex_buscado:
                            link.add(cara_set)
        return link
    
    def calcular_caracteristica_euler(self):
        """
        Calcula la característica de Euler de un complejo simplicial.

        La característica de Euler es un invariante topológico, una propiedad que permanece constante 
        a pesar de la deformación continua del espacio (como estirar o comprimir, pero sin rasgar o pegar). 
        Para un complejo simplicial K, la característica de Euler χ(K) se define como:
        
        χ(K) = Σ(-1)^k * s_k,
        
        donde s_k es el número de símplices de dimensión k en K. Esta suma alterna los signos entre las dimensiones.

        :return: La característica de Euler del complejo simplicial. Un número entero que es un invariante topológico del complejo.
        """
        num_simplices_por_dimension = {}
        for simplex, _ in self.simplices:
            # Calcula la dimensión de cada símplice, que es una menos que el número de sus vértices.
            dim = len(simplex) - 1
            # Cuenta el número de símplices en cada dimensión.
            if dim in num_simplices_por_dimension:
                num_simplices_por_dimension[dim] += 1
            else:
                num_simplices_por_dimension[dim] = 1

        # Realiza la suma alternante para calcular la característica de Euler.
        # La suma alterna los signos para cada dimensión (positivo para dimensiones pares, negativo para impares).
        caracteristica_euler = sum((-1)**dim * count for dim, count in num_simplices_por_dimension.items())
        return caracteristica_euler
        # Devuelve la característica de Euler, que es un indicador importante de la topología del complejo simplicial.

    
    def calcular_componentes_conexas(self):
        """
        Determina el número de componentes conexas en el espacio subyacente de un complejo simplicial.

        Una "componente conexa" en topología es una subsección de un espacio topológico donde cualquier 
        par de puntos dentro de esta subsección puede conectarse mediante un camino que permanezca 
        completamente dentro de la subsección. En el contexto de un complejo simplicial, una componente 
        conexa es un subconjunto de símplices donde cada par de símplices está conectado directa o 
        indirectamente a través de otros símplices del subconjunto.

        :return: El número de componentes conexas en el complejo simplicial. Este es un número entero que 
        representa cuántas "piezas" separadas e independientes forman el complejo.
        """
        visitados = set()  # Conjunto para marcar los vértices ya visitados en la búsqueda.
        num_componentes_conexas = 0  # Contador para el número de componentes conexas.

        def dfs(v):
            """
            Realiza una búsqueda en profundidad (DFS) para marcar todos los vértices conectados a 'v'.

            La búsqueda en profundidad es un algoritmo de recorrido de grafos que comienza en un vértice 
            y explora lo más profundamente posible a lo largo de cada rama antes de retroceder, lo que 
            es ideal para identificar componentes conexas.
            """
            if v in visitados:
                return  # Si el vértice ya ha sido visitado, no hay necesidad de explorarlo nuevamente.
            visitados.add(v)  # Marcar el vértice actual como visitado.
            for vecino in self.grafo[v]:
                dfs(vecino)  # Explorar recursivamente todos los vecinos del vértice actual.

        # Iniciar DFS en cada vértice no visitado para identificar todas las componentes conexas.
        for vertice in self.grafo:
            if vertice not in visitados:
                dfs(vertice)  # Iniciar DFS en este vértice.
                num_componentes_conexas += 1  # Cada llamada a DFS que inicia una nueva exploración indica una nueva componente conexa.

        return num_componentes_conexas
        # Devuelve el número total de componentes conexas encontradas en el complejo simplicial.

    
    def subcomplejo_por_filtracion(self, valor_filtracion):
        """
        Genera un subcomplejo simplicial basado en un valor de filtración dado.

        La "filtración" en topología computacional es un proceso que construye complejos simpliciales 
        incrementando estructuras sobre la base de algún criterio, como una medida de similitud, distancia, 
        densidad, etc. Este proceso permite analizar la estructura topológica del complejo a diferentes 
        'escalas' o niveles de detalle.

        :param valor_filtracion: El valor de filtración máximo que se utiliza como umbral para incluir 
        un símplice en el subcomplejo. Símplices con un valor de filtración más alto que este umbral 
        no se incluyen en el subcomplejo.

        :return: Un subcomplejo simplicial que contiene solo los símplices cuyo valor de filtración es 
        menor o igual que el valor especificado. Este subcomplejo representa una 'vista' del complejo 
        original a un nivel particular de filtración.
        """
        subcomplejo = ComplejoSimplicial()  # Inicializa un nuevo complejo simplicial vacío.

        for simplex, filtracion in self.simplices:
            # Itera sobre todos los símplices en el complejo simplicial junto con sus valores de filtración.
            if filtracion <= valor_filtracion:
                # Incluye el símplice en el subcomplejo si su valor de filtración es menor o igual al umbral dado.
                subcomplejo.añadir_simplex(simplex, filtracion)
                # Añade el símplice al subcomplejo.

        return subcomplejo
        # Devuelve el subcomplejo simplicial resultante, que es una subestructura del complejo original 
        # filtrada según el valor de filtración especificado.

    def calcular_matriz_borde(self, dimension):
        """
        Calcula la matriz borde de una dimensión dada en un complejo simplicial.

        En topología algebraica, una "matriz borde" es una representación matricial que describe cómo 
        los símplices de una dimensión específica se 'conectan' o 'limitan' con los símplices de la dimensión 
        inmediatamente inferior. Esta matriz es fundamental para calcular los grupos de homología de un complejo.

        :param dimension: La dimensión de los símplices para los cuales se calculará la matriz borde.

        :return: Una matriz (en forma de array NumPy) donde cada fila corresponde a un símplice de dimensión (dimension-1)
        y cada columna corresponde a un símplice de dimensión 'dimension'. Los elementos de la matriz son 0 o 1,
        indicando si el símplice de dimensión inferior es o no una cara del símplice de dimensión superior.
        """

        # Crear listas de símplices por dimensión
        simplices_d = list(self.caras_por_dimension(dimension))
        simplices_d_minus_1 = list(self.caras_por_dimension(dimension-1))

        # Crear la matriz borde inicializando todos los elementos a 0
        matriz_borde = np.zeros((len(simplices_d_minus_1), len(simplices_d)), dtype=int)

        # Rellenar la matriz borde
        for j, simplex in enumerate(simplices_d):
            for i, sub_simplex in enumerate(simplices_d_minus_1):
                # Si un simplex de dimensión inferior es una cara del simplex de dimensión superior,
                # coloca un 1 en la matriz en la posición correspondiente.
                if sub_simplex.issubset(simplex):
                    matriz_borde[i, j] = 1

        return matriz_borde
        # La matriz resultante es la matriz borde para la dimensión especificada.

    
    
    def calcular_matriz_borde_generalizada(self):
        """
        Calcula una matriz borde generalizada para un complejo simplicial filtrado.

        En un complejo simplicial filtrado, cada símplice se asocia con un valor de filtración, 
        que generalmente representa algún tipo de medida o umbral (como la distancia, el tiempo, 
        la densidad, etc.). Esta matriz borde generalizada tiene en cuenta no solo las relaciones 
        topológicas de los símplices, sino también cómo estas relaciones se desarrollan a través 
        de la filtración.

        :return: Una matriz (en realidad, una lista de listas debido a la naturaleza irregular 
        de los datos) que representa las relaciones borde entre los símplices en diferentes 
        dimensiones, considerando sus valores de filtración.
        """

        # Ordenar símplices por filtración y dimensión
        self.simplices.sort(key=lambda x: (x[1], len(x[0])))
        # Esto organiza los símplices primero por su valor de filtración y luego por su dimensión.

        # Crear listas de símplices organizadas por dimensión y filtración
        simplices_organizados = defaultdict(list)
        for simplex, filtracion in self.simplices:
            dim = len(simplex) - 1
            simplices_organizados[dim].append((simplex, filtracion))
            # Agrupa los símplices por su dimensión, manteniendo su valor de filtración.

        # Determinar el número máximo de dimensiones en el complejo
        max_dim = max(simplices_organizados.keys())

        # Crear la matriz borde generalizada
        bordes = []
        for dim in range(1, max_dim + 1):
            fila = []
            for simplex_d, filtracion_d in simplices_organizados[dim]:
                columna = []
                for simplex_d_minus_1, filtracion_d_minus_1 in simplices_organizados[dim - 1]:
                    # Verifica si un simplex de dimensión inferior es una cara del simplex de dimensión superior
                    # y si su valor de filtración es menor o igual al del simplex de dimensión superior.
                    if simplex_d_minus_1.issubset(simplex_d) and filtracion_d_minus_1 <= filtracion_d:
                        columna.append(1)
                    else:
                        columna.append(0)
                fila.append(columna)
            bordes.append(fila)

        return bordes
        # Devuelve la matriz borde generalizada, que refleja la estructura y evolución del complejo simplicial filtrado.


    def calcular_numeros_betti(self):
        """
        Calcula los números de Betti para el complejo simplicial.

        Los números de Betti son una secuencia de números enteros (β₀, β₁, β₂, ...) que representan 
        invariantes topológicos de un espacio. En particular, β₀ representa el número de componentes 
        conexas, β₁ el número de "agujeros" o ciclos de una dimensión, β₂ el número de "cavidades" 
        tridimensionales, y así sucesivamente.

        Returns:
        list: Los números de Betti (β₀, β₁, β₂, ...) del complejo simplicial.
        """

        st = gd.SimplexTree()  # Crear un árbol de símplices, una estructura de datos para almacenar el complejo simplicial.

        # Añadir todos los símplices al complejo simplicial
        for simplex, _ in self.simplices:
            st.insert(list(simplex))
            # Cada símplice se añade al árbol de símplices. Se asume que 'self.simplices' es una lista de tuplas,
            # donde cada tupla contiene un conjunto que representa un símplice y su valor de filtración.

        # Calcular la persistencia para todas las dimensiones
        st.persistence()
        # Este método calcula la homología persistente del complejo, que es necesaria para determinar los números de Betti.

        # Calcular y devolver los números de Betti
        return st.betti_numbers()
        # Devuelve una lista de los números de Betti. Cada elemento de la lista corresponde a un número de Betti para
        # una dimensión específica (empezando por la dimensión 0).



# Ejemplos de los primeros metodos de la clase
complejo = ComplejoSimplicial()

complejo.añadir_simplex({1, 2, 3}, 0.4)
complejo.añadir_simplex({2, 3, 4}, 0.45)
complejo.añadir_simplex({1, 2}, 0.6) # Añadiendo simplices

# Convertir frozenset a listas para una mejor visualización
todas_las_caras = [list(cara) for cara in complejo.obtener_todas_las_caras()]
caras_dimension_1 = [list(cara) for cara in complejo.caras_por_dimension(0)]

print("Todas las caras del complejo:", len(todas_las_caras))
print("\nCaras de dimensión 1:", caras_dimension_1)

# Calculando la estrella del símplice {1, 2}
estrella_simplex = [list(s) for s in complejo.calcular_estrella({2, 3})]
print("\nEstrella del símplice {1, 2}:", estrella_simplex)

# Calculando el link del símplice {2, 3}
link_simplex = [list(s) for s in complejo.calcular_link({2, 3})]

print('\nEl link del simplice {2,3} es:', link_simplex)


# Calcular la característica de Euler del complejo
caracteristica_euler = complejo.calcular_caracteristica_euler()
print('\nla caracteristica de euler para el complejo simplicial es:', caracteristica_euler)

# Calcular subcomplejo por filtracion con un valor dado
valor_filtracion = 0.5
subcomplejo = complejo.subcomplejo_por_filtracion(valor_filtracion)
simplices_en_subcomplejo = [simplex for simplex, _ in subcomplejo.simplices]
print('\nLos simplices del subcomplejo para un valor de filtración de 0.5 es:', simplices_en_subcomplejo)

# Ejemplo para matriz de borde
complejo2 = ComplejoSimplicial()
complejo2.añadir_simplex({0, 1}, 0.1)
complejo2.añadir_simplex({1, 2}, 0.1)
complejo2.añadir_simplex({2, 3}, 0.1)
complejo2.añadir_simplex({0, 3}, 0.1)
complejo2.añadir_simplex({0, 2}, 0.1)
complejo2.añadir_simplex({1, 3}, 0.1)
complejo2.añadir_simplex({0, 1, 2}, 0.2)
complejo2.añadir_simplex({0, 2, 3}, 0.2)

matriz_borde_0 = complejo2.calcular_matriz_borde(1)
matriz_borde_1 = complejo2.calcular_matriz_borde(2)

print("\nMatriz borde de dimensión 0 a 1:\n", matriz_borde_0)
print("\nMatriz borde de dimensión 1 a 2:\n", matriz_borde_1)

##################################################################
# Practica 2

class AlfaComplejo:
    def __init__(self, puntos):
        """
        Inicializa un AlfaComplejo con un conjunto dado de puntos.

        :param puntos: Una lista de coordenadas que representan los puntos en el espacio.
        """
        self.puntos = np.array(puntos)  # Convierte los puntos en un array de NumPy para facilitar el cálculo.
        self.alfa_complejo = []  # Inicializa una lista vacía para almacenar el complejo Alfa.

    def calcular_filtracion_alfa(self, r):
        """
        Calcula el complejo Alfa para un valor de radio 'r'.

        En un Alfa-complejo, para un valor específico de r, incluimos en el complejo simplicial aquellos 
        símplices cuyas esferas circunscritas tienen un radio menor o igual a r.

        :param r: El radio para calcular la filtración del Alfa-complejo.
        :return: Una lista de símplices (cada simplex representado como una lista de índices de puntos) 
                 que forman el Alfa-complejo para el valor dado de r.
        """
        if r < 0:
            # Un radio negativo no tiene sentido en el contexto de un Alfa-complejo, por lo que se devuelve una lista vacía.
            return []
        elif r == 0:
            # Cuando r es 0, el Alfa-complejo solo consiste en los puntos individuales.
            # Cada punto es considerado un simplex de 0 dimensiones.
            return [[i] for i in range(len(self.puntos))]
        else:
            # Para r positivo, se calcula el complejo simplicial de Delaunay, que es la base para formar el Alfa-complejo.
            try:
                delaunay = Delaunay(self.puntos)
                # Se obtienen los símplices de la triangulación de Delaunay, que forman la base del Alfa-complejo.
                # Cada simplex está representado como una lista de índices de puntos.
                return [list(simplice) for simplice in delaunay.simplices]
            except :
                # Manejo de error: si los puntos son colineales o coplanares, la triangulación de Delaunay no se puede calcular.
                # Esto puede ocurrir, por ejemplo, si todos los puntos están en una misma línea o en un mismo plano.
                print("Los puntos podrían ser colineales o coplanares, lo que impide calcular la triangulación de Delaunay.")
                return None
            
    def representar_graficamente(self, r, guardar_archivo=True, nombre_archivo="alfa_complejo.png"):
        """
        Representa gráficamente el alfa complejo para un valor de r dado y opcionalmente lo guarda.

        :param r: El valor de r para el cual se representa el alfa complejo.
        :param guardar_archivo: Un booleano que indica si se debe guardar la representación en el sistema de archivos.
        :param nombre_archivo: El nombre del archivo bajo el cual se guardará la representación gráfica.
        """
        filtracion = self.calcular_filtracion_alfa(r)
        
        if filtracion is None:
            print("No se puede representar gráficamente el alfa complejo debido a un error anterior.")
            return
        
        # Extraer los puntos del complejo simplicial
        puntos_x = self.puntos[:, 0]
        puntos_y = self.puntos[:, 1]
        
        # Configurar la gráfica
        plt.figure()
        plt.axis([min(puntos_x) - 1, max(puntos_x) + 1, min(puntos_y) - 1, max(puntos_y) + 1])
        plt.gca().set_aspect('equal')
        
        # Dibujar puntos
        plt.plot(puntos_x, puntos_y, 'o')
        
        # Dibujar aristas y triángulos
        if r > 0:
            for simplex in filtracion:
                simplex_puntos = self.puntos[list(simplex)]
                if len(simplex) == 2:
                    # Dibujar arista
                    plt.plot(simplex_puntos[:, 0], simplex_puntos[:, 1], 'k-')
                elif len(simplex) == 3:
                    # Dibujar triángulo
                    triangulacion = tri.Triangulation(simplex_puntos[:, 0], simplex_puntos[:, 1])
                    plt.triplot(triangulacion, 'k-')
        
        # Guardar o mostrar la figura
        if guardar_archivo:
            plt.savefig(nombre_archivo)
            print(f"La representación gráfica se ha guardado como {nombre_archivo}")
        else:
            plt.show()

        plt.close()

        
# Ejmeplo para los alfacomplejos coplanares:
puntos = [(1, 2), (3, 4), (5, 6), (2, 5)]  # Conjunto de puntos del plano
alfa_complejo = AlfaComplejo(puntos)

# Valores de r para calcular la filtración
valores_r = [-1, 0, 2]  # r negativo, r cero y r suficientemente grande

# Calcula la filtración de alfa complejos 
filtraciones = {r: alfa_complejo.calcular_filtracion_alfa(r) for r in valores_r}
filtraciones

print('las filtraciones para diferentes r para los puntos [(1, 2), (3, 4), (5, 6), (2, 5)] son :', filtraciones)

# Representación del alfa complejo 
alfa_complejo.representar_graficamente(2)



class Complejosimplicial2:
    def __init__(self, puntos):
        self.puntos = np.array(puntos)
    
    def calcular_filtracion_vietoris_rips(self, r):
        if r < 0:
            return []
        
        distancias = squareform(pdist(self.puntos))
        simplices = []
        
        for k in range(len(self.puntos) + 1):
            for simplex in combinations(range(len(self.puntos)), k):
                if simplex:  # Evita agregar el conjunto vacío
                    # Comprobamos que todos los pares de puntos en el simplex están a una distancia <= r
                    if all(distancias[i][j] <= r for i, j in combinations(simplex, 2)):
                        simplices.append(list(simplex))
        return simplices

# Ejemplo de uso:
puntos_ejemplo = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])  # Cuatro puntos formando un cuadrado
complejo = Complejosimplicial2(puntos_ejemplo)
filtracion_rips = complejo.calcular_filtracion_vietoris_rips(1.5)

# Las filtraciones de Vietoris-Rips se devuelven como listas de índices
print('la filtracion de vieto-rips para el conjunto de puntos [[0, 0], [1, 0], [0, 1], [1, 1]] es:', filtracion_rips)


##############################################################
# Práctica 3: Homología Simplicial


def forma_normal_smith_z2(matriz):
    M = np.array(matriz, dtype=np.int64) % 2  # Asegurar que la matriz esté en Z2
    filas, columnas = M.shape

    for i in range(min(filas, columnas)):
        # Encontrar un elemento no nulo en la columna i y moverlo a la posición (i, i)
        for j in range(i, filas):
            if M[j, i] == 1:
                # Intercambiar filas si es necesario
                if j != i:
                    M[[i, j]] = M[[j, i]]
                break
        else:
            continue  # Continuar si no se encontró ningún elemento no nulo
        
        # Hacer cero todos los elementos fuera de la diagonal con operaciones de fila y columna
        for j in range(filas):
            if j != i and M[j, i] == 1:
                M[j] = (M[j] + M[i]) % 2
        for k in range(columnas):
            if k != i and M[i, k] == 1:
                M[:, k] = (M[:, k] + M[:, i]) % 2
    
    return M

# Ejemplo de uso:
matriz_ejemplo = [
    [1, 1, 0],
    [0, 1, 1],
    [1, 0, 1]
]

forma_normal_smith = forma_normal_smith_z2(matriz_ejemplo)
print("La forma normal de Smith de la matriz es:\n", forma_normal_smith)

### Calculo de los numeros de Betti 

## Para el tetraedro
complejo_tetraedro = ComplejoSimplicial()
# Añadir símplices para formar un tetraedro completo
complejo_tetraedro.añadir_simplex({0, 1}, 0)
complejo_tetraedro.añadir_simplex({0, 2}, 0)
complejo_tetraedro.añadir_simplex({0, 3}, 0)
complejo_tetraedro.añadir_simplex({1, 2}, 0)
complejo_tetraedro.añadir_simplex({1, 3}, 0)
complejo_tetraedro.añadir_simplex({2, 3}, 0)
complejo_tetraedro.añadir_simplex({0, 1, 2}, 0)
complejo_tetraedro.añadir_simplex({0, 1, 3}, 0)
complejo_tetraedro.añadir_simplex({0, 2, 3}, 0)
complejo_tetraedro.añadir_simplex({1, 2, 3}, 0)

# Calcular los números de Betti
betti_numbers = complejo_tetraedro.calcular_numeros_betti()
print("\nNúmeros de Betti del tetraedro:", betti_numbers)

## Para el borde del tetraedro
complejo = ComplejoSimplicial()

# Añadir aristas para formar el borde del tetraedro
complejo.añadir_simplex({0, 1}, 0)
complejo.añadir_simplex({0, 2}, 0)
complejo.añadir_simplex({0, 3}, 0)
complejo.añadir_simplex({1, 2}, 0)
complejo.añadir_simplex({1, 3}, 0)
complejo.añadir_simplex({2, 3}, 0)

# Calcular los números de Betti para el borde del tetraedro
betti_numbers_borde_tetraedro = complejo.calcular_numeros_betti()
print('\nNúmeros de Betti para el borde del tetrahedro son:',betti_numbers_borde_tetraedro)

## Para el plano proyectivo
plano_proyectivo = ComplejoSimplicial()

# Añadir símplices para formar el plano proyectivo
plano_proyectivo.añadir_simplex({0, 1, 2}, 0)
plano_proyectivo.añadir_simplex({0, 1, 3}, 0)
plano_proyectivo.añadir_simplex({0, 2, 3}, 0)
plano_proyectivo.añadir_simplex({1, 2, 3}, 0)

# Calcular los números de Betti para el plano proyectivo
betti_numbers_plano_proyectivo = plano_proyectivo.calcular_numeros_betti()
print("\nNúmeros de Betti del plano proyectivo:", betti_numbers_plano_proyectivo)

## Para la botella de Klein
botella_klein = ComplejoSimplicial()

# Añadir símplices para formar la botella de Klein con filtración
botella_klein.añadir_simplex({0, 1, 2, 3}, 0)
botella_klein.añadir_simplex({0, 1, 2}, 1)
botella_klein.añadir_simplex({0, 1, 3}, 2)
botella_klein.añadir_simplex({0, 2, 3}, 3)
botella_klein.añadir_simplex({1, 2, 3}, 4)

# Calcular los números de Betti para la botella de Klein
betti_numbers_botella_klein = botella_klein.calcular_numeros_betti()
print("Números de Betti de la botella de Klein:", betti_numbers_botella_klein)


## Para el anillo
anillo = ComplejoSimplicial()

# Añadir símplices para formar el anillo con filtración
# Añadir símplices exteriores del anillo
anillo.añadir_simplex({0, 1}, 0.5)
anillo.añadir_simplex({1, 2}, 0.5)
anillo.añadir_simplex({2, 3}, 0.5)
anillo.añadir_simplex({3, 0}, 0.5)

# Añadir símplices interiores del anillo
anillo.añadir_simplex({0, 1, 2}, 0.3)
anillo.añadir_simplex({1, 2, 3}, 0.3)
anillo.añadir_simplex({0, 2, 3}, 0.3)
anillo.añadir_simplex({0, 1, 3}, 0.3)

# Calcular los números de Betti para el anillo
betti_numbers_anillo = anillo.calcular_numeros_betti()
print("\nNúmeros de Betti del anillo:", betti_numbers_anillo)

## Para el sombrero del asno
sombrero_asno = ComplejoSimplicial()

# Añadir símplices para formar el sombrero del asno con filtración
sombrero_asno.añadir_simplex({0, 1, 2, 3}, 0)
sombrero_asno.añadir_simplex({0, 1, 2}, 1)
sombrero_asno.añadir_simplex({0, 1, 3}, 2)
sombrero_asno.añadir_simplex({0, 2, 3}, 3)
sombrero_asno.añadir_simplex({1, 2, 3}, 4)

# Calcular los números de Betti para el sombrero del asno
betti_numbers_sombrero_asno = sombrero_asno.calcular_numeros_betti()
print("\nNúmeros de Betti del sombrero del asno:", betti_numbers_sombrero_asno)

## Para los alfa complejos
# Crear una instancia de ComplejoSimplicial para el primer alfa complejo
alfa_complejo1 = ComplejoSimplicial()

# Añadir símplices para el primer alfa complejo con filtración
alfa_complejo1.añadir_simplex({0, 1, 2}, 0)
alfa_complejo1.añadir_simplex({0, 1}, 1)
alfa_complejo1.añadir_simplex({1, 2}, 2)
alfa_complejo1.añadir_simplex({0, 2}, 3)

# Crear una instancia de ComplejoSimplicial para el segundo alfa complejo
alfa_complejo2 = ComplejoSimplicial()

# Añadir símplices para el segundo alfa complejo con filtración
alfa_complejo2.añadir_simplex({2, 3, 4}, 0)
alfa_complejo2.añadir_simplex({2, 3}, 1)
alfa_complejo2.añadir_simplex({3, 4}, 2)
alfa_complejo2.añadir_simplex({2, 4}, 3)

# Calcular los números de Betti para el primer alfa complejo
betti_numbers_alfa1 = alfa_complejo1.calcular_numeros_betti()
print("\nNúmeros de Betti del primer alfa complejo:", betti_numbers_alfa1)

# Calcular los números de Betti para el segundo alfa complejo
betti_numbers_alfa2 = alfa_complejo2.calcular_numeros_betti()
print("\nNúmeros de Betti del segundo alfa complejo:", betti_numbers_alfa2)

## Numeros de Betti con el algoritmo incremental
def calcular_numeros_betti(simplices):
    """
    Calcular los números de Betti b0 y b1 para un complejo simplicial en el plano utilizando el algoritmo incremental.

    Parámetros:
    simplices (lista de tuplas): Lista de símplices, donde cada símplice es una tupla de vértices.

    Devuelve:
    tupla: Una tupla que contiene los números de Betti b0 y b1.
    """
    # Inicialización de los números de Betti
    b0 = 0  # Número de componentes conectadas
    b1 = 0  # Número de agujeros unidimensionales (ciclos)

    # Funciones auxiliares
    def find_set(vertex, parent):
        """
        Encontrar el representante del conjunto que contiene el vértice dado.
        Implementa la parte 'find' del algoritmo Union-Find.
        """
        if parent[vertex] != vertex:
            parent[vertex] = find_set(parent[vertex], parent)
        return parent[vertex]

    def union_sets(u, v, parent, rank):
        """
        Fusionar los conjuntos que contienen u y v.
        Implementa la parte 'union' del algoritmo Union-Find.
        """
        u_root = find_set(u, parent)
        v_root = find_set(v, parent)

        if u_root != v_root:
            # Fusionar el conjunto más pequeño en el más grande para mantener el árbol poco profundo
            if rank[u_root] < rank[v_root]:
                parent[u_root] = v_root
            elif rank[u_root] > rank[v_root]:
                parent[v_root] = u_root
            else:
                parent[v_root] = u_root
                rank[u_root] += 1
            return True  # Indica que ocurrió una fusión
        return False

    # Extraer vértices de los símplices
    vertices = set()
    for simplex in simplices:
        vertices.update(simplex)
    
    # Inicializar las estructuras de Union-Find
    parent = {vertex: vertex for vertex in vertices}
    rank = {vertex: 0 for vertex in vertices}

    # Procesar cada símplice
    for simplex in simplices:
        if len(simplex) == 1:
            b0 += 1  # Contar cada vértice como una nueva componente conectada
        elif len(simplex) == 2:
            # Si la arista conecta dos componentes previamente separadas, disminuir b0
            if union_sets(simplex[0], simplex[1], parent, rank):
                b0 -= 1
            else:
                # La arista forma un ciclo, aumentar b1
                b1 += 1

    return b0, b1



# Definición de los complejos alpha

# Primer complejo alpha: Un cuadrado con una diagonal
# Vértices: 1, 2, 3, 4
# Aristas: (1,2), (2,3), (3,4), (4,1), (1,3)
alpha_complex_1 = [(1,), (2,), (3,), (4,), (1,2), (2,3), (3,4), (4,1), (1,3)]

# Segundo complejo alpha: Dos triángulos separados
# Vértices: 1, 2, 3, 4, 5, 6
# Aristas: (1,2), (2,3), (3,1), (4,5), (5,6), (6,4)
alpha_complex_2 = [(1,), (2,), (3,), (4,), (5,), (6,), 
                   (1,2), (2,3), (3,1), (4,5), (5,6), (6,4)]

# Calcular los números de Betti para ambos complejos
b0_1, b1_1 = calcular_numeros_betti(alpha_complex_1)
b0_2, b1_2 = calcular_numeros_betti(alpha_complex_2)

print('Los números de Betti para los alpha complejos por el metodo de incremental:',(b0_1, b1_1), (b0_2, b1_2))


###########################
## Practica 4

# Ejemplo de matriz de borde generalizada para un complejo simplicial

# Crear una instancia de ComplejoSimplicial
complejo_ = ComplejoSimplicial()

# Añadir algunos símplices con valores de filtración
complejo_.añadir_simplex({1}, 0.1)  # vértice
complejo_.añadir_simplex({2}, 0.2)  # vértice
complejo_.añadir_simplex({3}, 0.3)  # vértice
complejo_.añadir_simplex({1, 2}, 0.4)  # arista
complejo_.añadir_simplex({2, 3}, 0.5)  # arista
complejo_.añadir_simplex({1, 3}, 0.6)  # arista
complejo_.añadir_simplex({1, 2, 3}, 0.7)  # triángulo

# Calcular y mostrar la matriz borde generalizada
matriz_borde_generalizada = complejo_.calcular_matriz_borde_generalizada()
for fila in matriz_borde_generalizada:
    for columna in fila:
        print(columna)
    print("\n")

# Función para calular el low de la columna de una matriz

def calcular_low_columna(matriz, columna_index):
    """
    Calcula el 'low' de una columna específica en una matriz de NumPy.

    :param matriz: Matriz de NumPy.
    :param columna_index: Índice de la columna para calcular el 'low'.
    :return: Índice de la fila más baja con un valor no nulo, o -1 si la columna es cero.
    """
    # Revisar si la columna es completamente cero
    if np.all(matriz[:, columna_index] == 0):
        return -1

    # Encontrar el índice de la fila más baja con un valor no nulo
    filas_no_cero = np.where(matriz[:, columna_index] != 0)[0]
    return max(filas_no_cero)

# Ejemplo de uso:
matriz_ejemplo = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 1]])
columna_index = 2  # Para la tercera colummna
low = calcular_low_columna(matriz_ejemplo, columna_index)
print('Para la matriz\n',matriz_ejemplo)
print(f"El 'low' de la columna {columna_index} es: {low}")


## Reducción por columnas una matriz cuadrada según el algoritmo matricula de cálculo de persistencia.
def reducir_matriz(matriz):
    """
    Reduce una matriz cuadrada según el algoritmo de cálculo de persistencia.

    :param matriz: Matriz de NumPy a reducir.
    :return: Matriz reducida.
    """
    num_columnas = matriz.shape[1]
    lows = [-1] * num_columnas  # Inicializa todos los 'lows' como -1

    for j in range(num_columnas):
        while True:
            low_j = calcular_low_columna(matriz, j)
            if low_j == -1:
                break  # Si la columna es cero, no hay nada que reducir

            if lows[low_j] == -1:
                lows[low_j] = j  # Actualiza el 'low' de la fila
                break
            else:
                q = lows[low_j]
                matriz[:, j] = matriz[:, j] ^ matriz[:, q]  # Operación XOR para la reducción

    return matriz

# Ejemplo de uso
matriz_ejemplo = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 1]], dtype=int)
matriz_reducida = reducir_matriz(matriz_ejemplo)
print(f'\nLa matriz del ejemplo reducida es:\n{matriz_reducida}')

# Definición de la función calcular_diagrama_persistencia
def calcular_diagrama_persistencia(puntos, max_r, paso_r):
    # Crear un objeto AlfaComplejo a partir de los puntos de entrada
    alfa_complejo = AlfaComplejo(puntos)
    
    # Inicializar una lista para almacenar las filtraciones
    filtraciones = []
    
    # Generar las filtraciones en intervalos de r desde 0 hasta max_r con paso_r
    for r in np.arange(0, max_r, paso_r):
        filtraciones.append((r, alfa_complejo.calcular_filtracion_alfa(r)))
    
    # Crear un objeto SimplexTree para representar el complejo simplicial filtrado
    st = gd.SimplexTree()
    
    # Insertar los símplices en el SimplexTree con sus respectivas filtraciones
    for r, filtracion in filtraciones:
        for simplex in filtracion:
            st.insert(simplex, r)
    
    # Calcular y devolver el diagrama de persistencia del complejo filtrado
    return st.persistence()

# Definición de la función dibujar_diagrama_persistencia
def dibujar_diagrama_persistencia(diagrama, ax):
    # Extraer los nacimientos y muertes del diagrama de persistencia
    nacimientos = [p[1][0] for p in diagrama if p[1][1] != float('inf')]
    muertes = [p[1][1] for p in diagrama if p[1][1] != float('inf')]
    
    # Crear un gráfico de dispersión para visualizar los puntos de nacimiento y muerte
    ax.scatter(nacimientos, muertes)
    
    # Dibujar la línea diagonal punteada que representa la vida infinita
    ax.plot([0, max(muertes)], [0, max(muertes)], 'k--')
    
    # Configurar etiquetas y título del gráfico
    ax.set_xlabel('Nacimiento')
    ax.set_ylabel('Muerte')
    ax.set_title('Diagrama de Persistencia')

# Definición de la función dibujar_codigos_de_barras
def dibujar_codigos_de_barras(diagrama, ax):
    # Obtener el número de características en el diagrama
    numero_de_caracteristicas = len(diagrama)
    
    # Recorrer el diagrama y dibujar los códigos de barras para cada característica
    for i, (dim, (nacimiento, muerte)) in enumerate(diagrama):
        if muerte == float('inf'):
            # Si la muerte es infinita, establecerla como el máximo valor no infinito en el diagrama
            muerte = max([n for d, (n, m) in diagrama if m != float('inf')])
        
        # Dibujar un código de barras para la característica actual
        ax.plot([nacimiento, muerte], [i, i], 'b')
    
    # Configurar etiquetas y título del gráfico
    ax.set_xlabel('Valor de filtración')
    ax.set_ylabel('Características')
    ax.set_title('Códigos de Barras')
    
    # Ajustar las etiquetas en el eje y para que estén correctamente espaciadas
    ax.set_yticks(np.linspace(0, numero_de_caracteristicas - 1, 7))
    ax.set_ylim(-1, numero_de_caracteristicas)

# Generar puntos para diferentes curvas (circunferencia, figura ocho, elipse)
theta = np.linspace(0, 2*np.pi, 100)

# Circunferencia
circunferencia = np.array([[np.cos(t), np.sin(t)] for t in theta])

# Figura ocho
figura_ocho = np.array([[np.sin(t), np.sin(t) * np.cos(t)] for t in theta])

# Elipse
elipse = np.array([[2*np.cos(t), np.sin(t)] for t in theta])

# Función para visualizar cada curva, su diagrama de persistencia y códigos de barras
def visualizar_curva_y_analisis(curva, titulo):
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    # Dibujar la curva
    axs[0].plot(curva[:, 0], curva[:, 1])
    axs[0].set_title(titulo)
    axs[0].set_aspect('equal')

    # Calcular y dibujar diagrama de persistencia y códigos de barras
    diagrama = calcular_diagrama_persistencia(curva, max_r=2.0, paso_r=0.1)
    dibujar_diagrama_persistencia(diagrama, axs[1])
    dibujar_codigos_de_barras(diagrama, axs[2])

    plt.tight_layout()
    plt.savefig(f'{titulo}.png')

# Visualizar cada curva con su análisis
visualizar_curva_y_analisis(circunferencia, "Circunferencia")
visualizar_curva_y_analisis(figura_ocho, "Figura Ocho")
visualizar_curva_y_analisis(elipse, "Elipse")