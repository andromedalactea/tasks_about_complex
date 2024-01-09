from itertools import combinations
from collections import defaultdict

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


##################################################################
# Practica 2

from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

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