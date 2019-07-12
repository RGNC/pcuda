readme.txt: Este fichero tiene como objetivo explicar la estructura del algoritmo y el código.

El objetivo de este simulador es el de simular un P sistema usando membranas activas. La entrada
del algoritmo viene dada por un fichero binario que tiene codificado la configuración inicial
del sistema, y sus reglas.
Esta versión inicial utiliza una versión del algoritmo secuencial (ideado inicialmente por Ignacio
Pérez-Hurtado en su simulador en java de p-lingua), mediante dos fases:
1. Selección de reglas a ser ejecutadas.
2. Ejecución de las reglas seleccionadas en 1.

El objetivo de cada fichero es el siguiente:

- main.cpp: Contiene la función main, y es el encargado de leer los parámetros de entrada al
simulador (ver el código o escribir -help para saber cuales son), llamar al parser y llamar al
algoritmo secuencial.

- parserbin.h, parserbin.cpp: Contiene dos clases bien diferenciadas:
  * Parserbin: Lee, decodifica y devuelve la información contenida en un fichero binario de entrada.
  * Configuration: Almacena la configuración de un P sistema (administra el conjunto de reglas, de
  membranas, alfabeto, etiquetas, reglas seleccionadas, identificador de las membranas, etc). Al
  leer un fichero de entrada, Parserbin devuelve un objeto Configuration con toda la información.

- binvalues.h: Es un fichero donde se definen macros que configuran aspectos como la longitud de
cada campo de un fichero binario, etc.

- membrane.h: Define la estructura de la que se compone una Membrana (uso de la estructura primer hijo
y siguiente hermano, en vez de llevar un listado por cada membrana de todas sus hijas).

- multiset.h, multiset.cpp: Implementa un multiconjunto mediante una lista enlazada de los objetos que contiene (si
un objeto tiene multiplicidad 0, no está en la lista). Contiene métodos útiles para su uso posterior,
y uno para devolver, en forma de array de multiplicidades de todos los objetos, el contenido del
multiconjunto.

- ruleset.h, ruleset.cpp: Define clases útiles para almacenar y recorrer reglas:
	* Ruleset: Es un almacen de reglas. Guarda las reglas de forma ordenada: Diferencia primero por
	etiqueta, luego por carga, y por último, una lista de reglas ordenadas por tipo (primero las de
	evolución, y después la demás).
	* Rulelist: Es una clase iteradora de reglas. En vez de devolver directamente la lista de reglas
	asociada a cada etiqueta y carga, se devuelve esta clase que contiene métodos para iniciar una
	iteración, pasar al siguiente elemento, etc.
	* Selectedrulelist: Define la lista de reglas seleccionadas que será utilizada en la fase de selección.
	En cada nodo se almacena la regla seleccionada, la membrana que ha sido seleccionada y el número de veces
	que se selecciona tal regla. Además, contiene métodos para recorrer la lista iterativamente.
	
- sequential.h, sequential.cpp: Define las funciones que implementan el algoritmo del simulador (sequential_solve)
y sus dos fases: selección (select_rules) y ejecución (execute_rules).

- timestat/timestat.h: Una librería, creada y podificada por Miguel A. Martínez-del-Amor, para la obtención e
impresión del tiempo transcurrido dentro del programa.

/*$Id: readme.txt 2009-01-11 15:14:44 mdelamor $*/
