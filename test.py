import gudhi
import numpy as np

# Definir los vértices del borde del tetraedro como puntos en el espacio 3D
tetrahedron_vertices = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
]

# Crear un objeto CubicalComplex especificando los vértices
cubical_complex = gudhi.CubicalComplex(vertices=tetrahedron_vertices)

# Calcular la persistencia
cubical_complex.compute_persistence()

# Especificar los valores de nacimiento y muerte para los números persistentes de Betti
from_value = 0.0  # Valor de nacimiento
to_value = 1.0    # Valor de muerte (considerando el borde completo)

# Calcular los números de Betti persistentes del borde del tetraedro
persistent_betti_numbers = cubical_complex.persistent_betti_numbers(from_value, to_value)

# Imprimir los números de Betti persistentes
print("Números de Betti persistentes del borde del tetraedro:", persistent_betti_numbers)
