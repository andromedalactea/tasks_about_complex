¿Por qué se usa frozenset para el desarrollo?
Inmutabilidad:
Un frozenset es una versión inmutable de un conjunto (set) en Python. Una vez creado, no puede ser modificado (es decir, no puedes agregar o eliminar elementos).
Esta inmutabilidad es útil cuando se trabaja con estructuras de datos que no deben cambiar una vez definidas, como es el caso de los símplices en un complejo simplicial. Asegura la integridad de los datos a lo largo del tiempo.
Usabilidad como Elemento de un Conjunto:

En Python, los elementos de un conjunto (set) deben ser inmutables. Dado que los conjuntos (set) mismos son mutables, no pueden ser usados como elementos de otros conjuntos.
Sin embargo, frozenset, al ser inmutable, puede ser usado como un elemento dentro de otro conjunto. Esto es crucial cuando se desea tener un conjunto de conjuntos, como en el caso de almacenar las caras de un símplice, que son en sí mismas conjuntos de vértices.
Evitar Duplicados:

Al trabajar con complejos simpliciales, es importante asegurarse de que no haya duplicados en las caras. Los conjuntos (set y frozenset) son útiles porque automáticamente evitan duplicados.
Al convertir las caras a frozenset antes de almacenarlas en un conjunto, garantizamos que cada cara sea única en la colección.

Explicación del Método calcular_estrella:
Argumento simplex_buscado: Este es el símplice del cual queremos encontrar la estrella.
Bucle sobre self.simplices: El método recorre todos los símplices en el complejo.
Condición if simplex_buscado.issubset(simplex): Para cada símplice en el complejo, esta condición verifica si simplex_buscado es una cara de ese símplice. Si lo es, el símplice se añade al conjunto de la estrella.

¿Por que se usa defaultdict?
En el caso del método añadir_simplex, lo usamos para crear un conjunto asociado a cada vértice en el grafo del complejo simplicial. Cuando agregamos una arista, que es un simplex de dimensión 1 representado por un par de vértices (u, v), queremos asegurarnos de que ambos vértices u y v estén en el grafo y tengan un conjunto asociado con ellos para mantener el registro de sus vértices adyacentes. Al usar defaultdict(set), si intentamos acceder al conjunto de adyacencias de un vértice que no ha sido visto antes, defaultdict automáticamente crea un nuevo conjunto vacío como valor por defecto, lo que evita que tengamos que comprobar manualmente si el vértice existe en el diccionario antes de agregar su vecino

Este comportamiento hace que el código sea más limpio y menos propenso a errores, ya que no necesitamos insertar condicionales para inicializar las claves la primera vez que las encontramos. En resumen, defaultdict simplifica la creación y la gestión del grafo de adyacencias en el complejo simplicial.

¿Por que usar la busqueda DFS?
La búsqueda en profundidad (DFS, por sus siglas en inglés "Depth-First Search") es un algoritmo para recorrer o buscar en un árbol, árbol de búsqueda o grafo. Este método comienza en la raíz (seleccionando algún vértice arbitrario como raíz en el caso de un grafo) y explora tanto como sea posible a lo largo de cada rama antes de retroceder.

La implementación de DFS en el contexto de un complejo simplicial es útil por las siguientes razones:

Identificación de Componentes Conexas: En un grafo, una componente conexa es un subgrafo en el cual cualquier dos vértices están conectados entre sí por caminos, y al cual no están conectados vértices de ningún otro subgrafo. DFS puede ser utilizado para identificar estas componentes conexas explorando el grafo desde un vértice inicial y marcando todos los vértices alcanzables antes de pasar a un nuevo conjunto de vértices no marcados y repetir el proceso.

Eficiencia: DFS es eficiente en términos de memoria porque su almacenamiento es proporcional al máximo tamaño de la pila de llamadas recursivas, que en el peor de los casos es igual al número de vértices en el grafo.

Simplicidad: El algoritmo es conceptualmente simple y se puede implementar de manera recursiva o iterativa.

En el contexto de la topología y los complejos simpliciales, usar DFS para calcular el número de componentes conexas permite determinar cuántas partes distintas del complejo no están conectadas entre sí por símplices de dimensión uno (aristas). Por ejemplo, si tienes varios triángulos en un espacio, algunos de los cuales comparten aristas, DFS te ayudará a determinar cuántos grupos de triángulos conectados hay. Esto es vital para entender la estructura subyacente y las propiedades de conectividad del espacio que el complejo simplicial está modelando.
