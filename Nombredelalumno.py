from itertools import combinations
from collections import defaultdict
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import gudhi as gd 
# Parte 1 y 3
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
        Calcula el conjunto de todas las caras del complejo simplicial.

        :return: Un conjunto con todas las caras del complejo.
        """
        todas_las_caras = set()
        for simplex, _ in self.simplices:
            for i in range(len(simplex)):
                for cara in combinations(simplex, i + 1):
                    todas_las_caras.add(frozenset(cara))
        return todas_las_caras

    def caras_por_dimension(self, dimension):
        """
        Calcula el conjunto de todas las caras de una dimensión dada.

        :param dimension: La dimensión de las caras a calcular.
        :return: Un conjunto con todas las caras de la dimensión especificada.
        """
        return {cara for cara in self.obtener_todas_las_caras() if len(cara) - 1 == dimension}

    def calcular_estrella(self, simplex_buscado):
        """
        Calcula la estrella de un símplice dado.

        :param simplex_buscado: El símplice para el cual se calculará la estrella.
        :return: Un conjunto con todos los símplices que forman la estrella del símplice dado.
        """
        simplex_buscado = frozenset(simplex_buscado)
        estrella = set()
        for simplex, _ in self.simplices:
            if simplex_buscado.issubset(simplex):
                estrella.add(frozenset(simplex))
        return estrella
    
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

        La característica de Euler χ(K) se define como la suma alternante del número
        de símplices en cada dimensión, es decir, χ(K) = Σ(-1)^k * s_k,
        donde s_k es el número de símplices de dimensión k en K.

        :return: La característica de Euler del complejo simplicial.
        """
        num_simplices_por_dimension = {}
        for simplex, _ in self.simplices:
            # La dimensión de un simplex es una menos que el número de sus vértices
            dim = len(simplex) - 1
            if dim in num_simplices_por_dimension:
                num_simplices_por_dimension[dim] += 1
            else:
                num_simplices_por_dimension[dim] = 1

        # Calcular la suma alternante
        caracteristica_euler = sum((-1)**dim * count for dim, count in num_simplices_por_dimension.items())
        return caracteristica_euler
    
    def calcular_componentes_conexas(self):
        """
        Determina el número de componentes conexas del espacio subyacente al complejo simplicial.

        :return: El número de componentes conexas del complejo simplicial.
        """
        visitados = set()
        num_componentes_conexas = 0

        def dfs(v):
            """
            Realiza una búsqueda en profundidad (DFS) para marcar todos los vértices conectados a 'v'.
            """
            if v in visitados:
                return
            visitados.add(v)
            for vecino in self.grafo[v]:
                dfs(vecino)

        # Iniciar DFS en cada vértice no visitado
        for vertice in self.grafo:
            if vertice not in visitados:
                dfs(vertice)
                num_componentes_conexas += 1

        return num_componentes_conexas
    
    def subcomplejo_por_filtracion(self, valor_filtracion):
        """
        Calcula el subcomplejo simplicial formado por todos los símplices cuyo valor de filtración
        es menor o igual que un valor dado.

        :param valor_filtracion: El valor de filtración máximo para incluir un simplex en el subcomplejo.
        :return: Un nuevo complejo simplicial que contiene solo los símplices con filtración
                 menor o igual que valor_filtracion.
        """
        subcomplejo = ComplejoSimplicial()
        for simplex, filtracion in self.simplices:
            if filtracion <= valor_filtracion:
                subcomplejo.añadir_simplex(simplex, filtracion)
        return subcomplejo
    def calcular_matriz_borde(self, dimension):
        # Crear listas de símplices por dimensión
        simplices_d = list(self.caras_por_dimension(dimension))
        simplices_d_minus_1 = list(self.caras_por_dimension(dimension-1))

        # Crear la matriz borde
        matriz_borde = np.zeros((len(simplices_d_minus_1), len(simplices_d)), dtype=int)

        # Rellenar la matriz borde
        for j, simplex in enumerate(simplices_d):
            for i, sub_simplex in enumerate(simplices_d_minus_1):
                if sub_simplex.issubset(simplex):
                    matriz_borde[i, j] = 1

        return matriz_borde
    def calcular_matriz_borde_generalizada(self):
        """
        Calcula la matriz borde generalizada para un complejo simplicial filtrado.

        :return: Una matriz que representa la relación borde entre símplices en diferentes dimensiones, teniendo en cuenta la filtración.
        """
        # Ordenar símplices por filtración y dimensión
        self.simplices.sort(key=lambda x: (x[1], len(x[0])))

        # Crear listas de símplices organizadas por dimensión y filtración
        simplices_organizados = defaultdict(list)
        for simplex, filtracion in self.simplices:
            dim = len(simplex) - 1
            simplices_organizados[dim].append((simplex, filtracion))

        # Determinar el número máximo de dimensiones
        max_dim = max(simplices_organizados.keys())

        # Crear la matriz borde generalizada
        bordes = []
        for dim in range(1, max_dim + 1):
            fila = []
            for simplex_d, filtracion_d in simplices_organizados[dim]:
                columna = []
                for simplex_d_minus_1, filtracion_d_minus_1 in simplices_organizados[dim - 1]:
                    if simplex_d_minus_1.issubset(simplex_d) and filtracion_d_minus_1 <= filtracion_d:
                        columna.append(1)
                    else:
                        columna.append(0)
                fila.append(columna)
            bordes.append(fila)

        return bordes

    def calcular_numeros_betti(self):
        """
        Calcula los números de Betti para cada dimensión del complejo simplicial.

        :return: Lista con los números de Betti.
        """
        max_dim = self.calcular_dimension()
        betti_numbers = []

        # Inicializar rangos de las matrices bordes con 0
        rangos_bordes = [0]

        # Calcular rangos de las matrices bordes para todas las dimensiones
        for dim in range(max_dim + 1):
            matriz_borde = self.calcular_matriz_borde(dim)
            rango_borde = np.linalg.matrix_rank(matriz_borde) if matriz_borde.size > 0 else 0
            rangos_bordes.append(rango_borde)

        # Calcular los números de Betti usando los rangos de las matrices bordes
        for dim in range(max_dim + 1):
            num_simplices_dim_actual = len(self.caras_por_dimension(dim))
            betti_num = num_simplices_dim_actual - rangos_bordes[dim]
            if dim < max_dim:
                betti_num -= rangos_bordes[dim + 1]
            betti_numbers.append(betti_num)

        return betti_numbers

    # Métodos adicionales (obtener_caras, caras_por_dimension, estrella, link, caracteristica_euler,
    # numero_componentes_conexas, ordenamiento, subcomplejo_por_filtracion) se implementarían aquí.

# Ejemplo de uso
complejo = ComplejoSimplicial()
complejo.añadir_simplex({1, 2, 3}, 0.4)
complejo.añadir_simplex({2, 3, 4}, 0.45)
complejo.añadir_simplex({1, 2}, 0.6)

# Convertir frozenset a listas para una mejor visualización
todas_las_caras = [list(cara) for cara in complejo.obtener_todas_las_caras()]
caras_dimension_1 = [list(cara) for cara in complejo.caras_por_dimension(0)]

print("Todas las caras del complejo:", len(todas_las_caras))
print("Caras de dimensión 1:", caras_dimension_1)

# Calculando la estrella del símplice {1, 2}
estrella_simplex = [list(s) for s in complejo.calcular_estrella({2, 3})]
print("Estrella del símplice {1, 2}:", estrella_simplex)

# Calculando el link del símplice {2, 3}
link_simplex = [list(s) for s in complejo.calcular_link({2, 3})]
print('El link del simplice {2,3} es:', link_simplex)

# Calcular la característica de Euler del complejo
caracteristica_euler = complejo.calcular_caracteristica_euler()
print('la caracteristica de euler para el complejo simplicial es:', caracteristica_euler)

# Calcular subcomplejo por filtracion con un valor dado
valor_filtracion = 0.5
subcomplejo = complejo.subcomplejo_por_filtracion(valor_filtracion)
simplices_en_subcomplejo = [simplex for simplex, _ in subcomplejo.simplices]
print('Los simplices del subcomplejo para un valor de filtración de 0.5 es:', simplices_en_subcomplejo)

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

print("Matriz borde de dimensión 0 a 1:\n", matriz_borde_0)
print("Matriz borde de dimensión 1 a 2:\n", matriz_borde_1)

##################################################################
# Practica 2


class AlfaComplejo:
    def __init__(self, puntos):
        self.puntos = np.array(puntos)
        self.alfa_complejo = []

    def calcular_filtracion_alfa(self, r):
        if r < 0:
            return []
        elif r == 0:
            # Cada punto es un simplex de 0 dimensiones por sí mismo
            return [[i] for i in range(len(self.puntos))]
        else:
            # Calcular el complejo simplicial de Delaunay para r suficientemente grande
            try:
                delaunay = Delaunay(self.puntos)
                # Obtener los índices de los puntos que forman cada simplex del complejo de Delaunay
                return [list(simplice) for simplice in delaunay.simplices]
            except :
                # Manejo de error si los puntos son colineales o coplanares
                
                print("Los puntos podrían ser colineales o coplanares, lo que impide calcular la triangulación de Delaunay.")
                return None
            
    def representar_graficamente(self, r, guardar_archivo=False, nombre_archivo="alfa_complejo.png"):
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

        
# Ejemplo de uso:
puntos = [(1, 2), (3, 4), (5, 6), (2, 5)]  # Conjunto de puntos del plano
alfa_complejo = AlfaComplejo(puntos)

# Valores de r para calcular la filtración
valores_r = [-1, 0, 2]  # r negativo, r cero y r suficientemente grande

# Intentamos calcular la filtración de alfa complejos nuevamente
filtraciones = {r: alfa_complejo.calcular_filtracion_alfa(r) for r in valores_r}
filtraciones

print('las filtraciones para diferentes r para los puntos [(1, 2), (3, 4), (5, 6), (2, 5)] son :', filtraciones)

# Representación del alfa complejo 
alfa_complejo.representar_graficamente(2)

from scipy.spatial.distance import pdist, squareform

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

# Añadir vértices (0-simplices)
for v in range(1, 5):
    complejo_tetraedro.añadir_simplex({v}, filtracion=0)

# Añadir aristas (1-simplices)
aristas = [{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}]
for arista in aristas:
    complejo_tetraedro.añadir_simplex(arista, filtracion=1)

# Añadir caras (2-simplices)
caras = [{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}]
for cara in caras:
    complejo_tetraedro.añadir_simplex(cara, filtracion=2)

# Calcular y mostrar los números de Betti
numeros_betti_tetraedro = complejo_tetraedro.calcular_numeros_betti()
print("Números de Betti del tetraedro:", numeros_betti_tetraedro)

## Para el borde del tetraedro
complejo_borde_tetraedro = ComplejoSimplicial()

# Añadir vértices (0-simplices)
for v in range(1, 5):
    complejo_borde_tetraedro.añadir_simplex({v}, filtracion=0)

# Añadir aristas (1-simplices)
aristas = [{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}]
for arista in aristas:
    complejo_borde_tetraedro.añadir_simplex(arista, filtracion=1)

# Añadir caras (2-simplices)
caras = [{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}]
for cara in caras:
    complejo_borde_tetraedro.añadir_simplex(cara, filtracion=2)

# Calcular y mostrar los números de Betti
numeros_betti_borde_tetraedro = complejo_borde_tetraedro.calcular_numeros_betti()
print("Números de Betti del borde del tetraedro:", numeros_betti_borde_tetraedro)


## Para el plano proyectivo
complejo_plano_proyectivo = ComplejoSimplicial()

# Añadir vértices (0-simplices)
for v in range(1, 3):
    complejo_plano_proyectivo.añadir_simplex({v}, filtracion=0)

# Añadir aristas (1-simplices)
aristas = [{1, 2}]
for arista in aristas:
    complejo_plano_proyectivo.añadir_simplex(arista, filtracion=1)

# Añadir el 2-simplex (cara)
complejo_plano_proyectivo.añadir_simplex({1, 2}, filtracion=2)

# Calcular y mostrar los números de Betti
numeros_betti_plano_proyectivo = complejo_plano_proyectivo.calcular_numeros_betti()
print("Números de Betti del plano proyectivo:", numeros_betti_plano_proyectivo)


## Para la botella de Klein
complejo_botella_klein = ComplejoSimplicial()

# Añadir vértices (0-simplices)
for v in range(1, 5):
    complejo_botella_klein.añadir_simplex({v}, filtracion=0)

# Añadir aristas (1-simplices)
aristas = [{1, 2}, {2, 3}, {3, 4}, {4, 1}, {1, 3}, {2, 4}]
for arista in aristas:
    complejo_botella_klein.añadir_simplex(arista, filtracion=1)

# Añadir caras (2-simplices)
caras = [{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}]
for cara in caras:
    complejo_botella_klein.añadir_simplex(cara, filtracion=2)

# Añadir el 3-simplex (volumen)
complejo_botella_klein.añadir_simplex({1, 2, 3, 4}, filtracion=3)

# Calcular y mostrar los números de Betti
numeros_betti_botella_klein = complejo_botella_klein.calcular_numeros_betti()
print("Números de Betti de la botella de Klein:", numeros_betti_botella_klein)



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
columna_index = 2  # Por ejemplo, elige la tercera columna
low = calcular_low_columna(matriz_ejemplo, columna_index)
print(f"El 'low' de la columna {columna_index} es: {low}")

# Reducción por columnas una matriz cuadrada según el algoritmo matricula de cálculo de persistencia.


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
print(f'La matriz del ejemplo reducida es:\n{matriz_reducida}')

# Definición de la función calcular_diagrama_persistencia
def calcular_diagrama_persistencia(puntos, max_r, paso_r):
    alfa_complejo = AlfaComplejo(puntos)
    filtraciones = []
    for r in np.arange(0, max_r, paso_r):
        filtraciones.append((r, alfa_complejo.calcular_filtracion_alfa(r)))
    st = gd.SimplexTree()
    for r, filtracion in filtraciones:
        for simplex in filtracion:
            st.insert(simplex, r)
    return st.persistence()

# Definición de la función dibujar_diagrama_persistencia
def dibujar_diagrama_persistencia(diagrama, ax):
    nacimientos = [p[1][0] for p in diagrama if p[1][1] != float('inf')]
    muertes = [p[1][1] for p in diagrama if p[1][1] != float('inf')]
    ax.scatter(nacimientos, muertes)
    ax.plot([0, max(muertes)], [0, max(muertes)], 'k--')
    ax.set_xlabel('Nacimiento')
    ax.set_ylabel('Muerte')
    ax.set_title('Diagrama de Persistencia')

# Definición de la función dibujar_codigos_de_barras
def dibujar_codigos_de_barras(diagrama, ax):
    numero_de_caracteristicas = len(diagrama)
    for i, (dim, (nacimiento, muerte)) in enumerate(diagrama):
        if muerte == float('inf'):
            muerte = max([n for d, (n, m) in diagrama if m != float('inf')])
        ax.plot([nacimiento, muerte], [i, i], 'b')
    ax.set_xlabel('Valor de filtración')
    ax.set_ylabel('Características')
    ax.set_title('Códigos de Barras')
    # Ajustar etiquetas en el eje y para que haya 7 correctamente espaciadas
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
    plt.show()

# Visualizar cada curva con su análisis
visualizar_curva_y_analisis(circunferencia, "Circunferencia")
visualizar_curva_y_analisis(figura_ocho, "Figura Ocho")
visualizar_curva_y_analisis(elipse, "Elipse")