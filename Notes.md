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