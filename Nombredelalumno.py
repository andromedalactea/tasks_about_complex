from itertools import combinations

class ComplejoSimplicial:
    def __init__(self):
        """
        Inicializa un complejo simplicial vacío.
        """
        self.simplices = []  # Lista para almacenar los símplices y sus valores de filtración

    def añadir_simplex(self, simplex, filtracion):
        """
        Añade un simplex al complejo simplicial y lo ordena.

        :param simplex: Un conjunto que representa un simplex (e.g., {1, 2, 3}).
        :param filtracion: Valor flotante asociado al simplex.
        """
        self.simplices.append((simplex, filtracion))
        self.simplices.sort(key=lambda x: (x[1], len(x[0])))

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

    # Métodos adicionales (obtener_caras, caras_por_dimension, estrella, link, caracteristica_euler,
    # numero_componentes_conexas, ordenamiento, subcomplejo_por_filtracion) se implementarían aquí.

# Ejemplo de uso
complejo = ComplejoSimplicial()
complejo.añadir_simplex({1, 2, 3}, 0.5)
complejo.añadir_simplex({2, 3, 4}, 0.6)
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