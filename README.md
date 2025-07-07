# Universidad Nacional de Colombia

## Facultad de Ciencias

## Programa de Ciencias de la Computacion

# Analisis de Bifurcaciones

## Presentado por: **Jacobo Ayala Giraldo**

## Materia: **Modelado y Simulacion**

## Fenomeno

Se desarrolla un análisis matemático de la dinámica de redes arbitrariamente pequeñas 
compuestas por poblaciones homogéneas de neuronas excitatorias e inhibitorias de tasa de 
disparo. Se estudian las bifurcaciones locales de su actividad neuronal con un enfoque que 
es tratable analiticamente , y se determina de manera numérica las bifurcaciones globales. 
Se encuentra que para una inhibición fuerte, estas redes dan lugar a unas dinámicas muy 
complejas, que emergen a través de la ruptura espontánea de simetría, es decir: 
Estado inicial simétrico: Las neuronas en la red inicialmente pueden tener estados 
equivalentes o simétricos 
Ruptura de simetría: Bajo ciertas condiciones (como inhibición fuerte), el sistema "decide" 
que algunas neuronas se comportan diferente a otras, creando patrones asimétricos de 
actividad 
Consecuencias: Esto genera múltiples soluciones ramificantes - es decir, diferentes 
patrones de actividad neuronal que el sistema puede adoptar, cada uno rompiendo la 
simetría de manera distinta 
En el cerebro, las neuronas reciben constantemente ‘ruido’, señales aleatorias no 
relacionadas con la información que estoy procesando. Cuando las neuronas son muy 
similares (homogéneas), tienden a responder de manera parecida tanto a las señales útiles 
como al ruido. 
Si muchas neuronas responden igual al ruido, este se amplifica y puede ‘tapar’ la señal real 
que el cerebro está tratando de procesar. 
Es allí donde nace la heterogeneidad funcional regulada, es decir, que las neuronas se 
especializan y responden de manera diferente a los mismos estímulos. Cada neurona 
desarrolla su propio ‘perfil’ de respuesta. La ruptura espontánea de la simetría hace que 
esta especialización ocurra de manera automática.

Este mecanismo es especialmente importante en redes pequeñas porque: 
* Las técnicas de campo medio (para redes grandes) no capturan este efecto 
* En redes pequeñas, cada neurona "cuenta más", así que la especialización tiene 
mayor impacto 
* Es un mecanismo que emerge naturalmente del tamaño intermedio del sistema

A escala mesoscópica, el cerebro se describe a menudo como una colección de masas 
neuronales, es decir, poblaciones neuronales homogéneas dentro de una columna cortical 
[30].Esta clase de modelos puede estudiarse usando la teoría de campo medio [35]. Con 
esta teoría aproximamos el comportamiento de redes grandes, y es útil para estudiar la 
actividad masiva de unos pocos miles de neuronas, lo que constituye un límite superior de 
las descripciones mesoscópicas de circuitos neuronales. 
Los modelos matemáticos de redes mesoscópicas deben poder predecir patrones que se 
observen en estos registros LFP/MEG/EEG.  Sin embargo, la teoría de campo medio 
proporciona sólo una aproximación del comportamiento real del sistema, y por lo tanto 
puede descuidar fenómenos importantes, como diferencias cualitativas en las transiciones 
entre regímenes estáticos y el caos [39]. Para ir más allá del límite de campo medio, se han 
desarrollado varias técnicas matemáticas para cuantificar los efectos de tamaño finito. Estos 
métodos de tamaño finito (como la aproximación de ruido lineal [41], el enfoque funcional de 
densidad [42], la teoría de grandes desviaciones [43] y los métodos de integral de 
trayectoria [44]) típicamente sólo pueden aplicarse a circuitos mesoscópicos compuestos 
por un número finito pero grande de neuronas. 
Los métodos para el análisis de redes hechas de unas pocas decenas de neuronas, que 
representan el límite inferior de la escala mesoscópica, siguen siendo limitados, el propósito 
de este trabajo es hacer progreso en la metodología matemática para estudiar la dinámica 
de tales redes.  El análisis de la dinámica de redes neuronales pequeñas fue pionero por 
Beer, quien estudió las bifurcaciones de redes de tamaño arbitrario con suposiciones 
altamente simétricas sobre la fuerza de los pesos sinápticos [48], y a través de 
aproximaciones asintóticas de las variedades de bifurcación [4]. En este artículo se extiende 
su análisis a una red más biológicamente plausible de tamaño arbitrario derivando 
expresiones exactas de las variedades de bifurcación, y con restricciones menos rígidas 
sobre los pesos sinápticos. El modelo se caracteriza por pesos homogéneos y 
arbitrariamente fuertes entre las poblaciones. Luego, realizamos un análisis numérico de las 
bifurcaciones globales que emergen al variar las corrientes de entrada externas (es decir, 
los estímulos) a la red y la fuerza de la inhibición, e introducimos una teoría matemática que 
describe las bifurcaciones locales analíticamente. Encontramos diferencias cualitativas con 
la aproximación de campo medio cuando el sistema tiene pesos sinápticos inhibitorios 
fuertes. En este caso, a través de un fenómeno de ruptura espontánea de simetría, la red 
neuronal experimenta una bifurcación especial conocida como punto de bifurcación o punto 
de ramificación [49], del cual emergen múltiples soluciones de las ecuaciones neuronales. 
En las nuevas ramas, pueden ocurrir nuevas bifurcaciones, enriqueciendo 
considerablemente la complejidad del diagrama de bifurcaciones de la red neuronal. Esta 
dinámica no es revelada por la aproximación de campo medio. 
Las suposiciones que se hacen sobre el modelo de circuitos neuronales de tamaño finito, es 
estructura en la siguiente figura. 

<image src="image.png" alt="Modelo de circuitos neuronales">

Se realiza el estudio numerico y analitico considerando el caso de dos poblaciones 
neronales de neuronas excitatorias e inhibitorias respectivamente. Las poblaciones 
contienen un numero finito arbitrario de neuronas, que estan conectadas entre si a traves de 
conexiones sinapticas no simetricas con pesos arbitrariemente fuertes. Para hacer la red 
anailiticamente tratable, se hacen algunas suposiciones simplificadoras. Se asume que las 
neuronas en cada poblacion tienen parametros homogeneos, que las redes neuronales son 
fully-connected y que los retrasos axonales son despreciables.

## Analisis analitico

### Ecuaciones del modelo

El modelo de red neuronal se describe mediante las siguientes ecuaciones diferenciales:

```math
\frac{dV_{i}(t)}{dt} = -\frac{1}{\tau_{i}}V_{i}(t) + \frac{1}{M_{i}}\sum_{j=0}^{N-1} J_{ij}\phi(V_{j}(t)) + I_{i}, \quad i=0, ..., N-1 \quad (1)
```
Donde:
* $N \geq 4$ representa el número de neuronas en la red
* $V_i(t)$ es el potencial de membrana de la neurona $i$ en el tiempo $t$
* $\tau_i$ es la constante de tiempo de la neurona $i$
* $M_i$ es el número de conexiones entrantes a la neurona $i$
* $J_{ij}$ es el peso sináptico de la conexión de la neurona $j$ a la neurona $i$
* $\phi(V_j(t))$ es la función de respuesta neuronal de la neurona $j$
* $I_i$ es la corriente de entrada externa a la neurona $i$

Típicamente se requieren funciones de activación lineales por tramos para obtener 
resultados analíticos [11, 55, 56]. Sin embargo, aquí mostraremos cómo obtener 
expresiones explícitas para los puntos de equilibrio y las bifurcaciones locales de la red 
cuando A(V) es una función suave (es decir, infinitamente diferenciable). En más detalle, 
consideraremos funciones de activación en forma de S (es decir, sigmoidales), ya que son 
biológicamente realistas [6, 54]. Particularmente conveniente desde el punto de vista 
matemático es la llamada función de activación algebraica, que se define como sigue:

```math
\phi_{j}(V) = \frac{v_{j}^{max}}{2} \left[ 1 + \frac{\frac{\Lambda_{j}}{2}(V - V_{j}^{T})}{\sqrt{1 + \frac{\Lambda_{j}^{2}}{4}(V - V_{j}^{T})^{2}}} \right]
```

Donde:
* $v_{j}^{max}$ es el valor máximo de la función de activación para la neurona $j$
* $V_{j}^{T}$ es el umbral de activación de la neurona $j$ (potencial de membrana en el cual la neurona comienza a activarse)
* $\Lambda_{j}$ es el parámetro de pendiente que controla qué tan abrupta es la transición de la función de activación alrededor del umbral


Para hacer nuestro análisis analíticamente tratable, suponemos que todos los parámetros 
del sistema son homogéneos para cada población de neuronas considerada. Así, estos 
parámetros serán indexados solo a nivel de población. Definiendo $N_E$ ($N_I$) como el tamaño 
de la población excitatoria (inhibitoria) (con $N_E, N_I \geq 2$), e indexando las neuronas de la 
población excitatoria como $i = 0, ..., N_E - 1$ y las neuronas inhibitorias como $i = N_E, ..., N - 1$ 
(con $N = N_E + N_I$), la matriz de conectividad sináptica se estructura como sigue:

```math
J = 
\begin{bmatrix} 
\mathcal{J}_{EE} & \mathcal{J}_{EI} \\ 
\mathcal{J}_{IE} & \mathcal{J}_{II} 
\end{bmatrix}, 
\quad 
\mathcal{J}_{\alpha\beta} = 
\begin{cases} 
J_{\alpha\alpha}(\Pi_{N_{\alpha}} - \text{Id}_{N_{\alpha}}), & \text{for } \alpha = \beta \\ 
J_{\alpha\beta}\Pi_{N_{\alpha}, N_{\beta}}, & \text{for } \alpha \ne \beta 
\end{cases}
```

Donde:
* $J_{EE}$ es la matriz de conectividad dentro de la población excitatoria
* $J_{EI}$ es la matriz de conectividad desde la población inhibitoria hacia la excitatoria  
* $J_{IE}$ es la matriz de conectividad desde la población excitatoria hacia la inhibitoria
* $J_{II}$ es la matriz de conectividad dentro de la población inhibitoria

Claramente tenemos $J_{EE}, J_{IE} > 0$, y $J_{II}, J_{EI} < 0$, lo que también significa que $M_E = M_I = N - 1$

```math
\mathbf{I}_\alpha = I_\alpha \mathbf{1}_{N_\alpha}
```

Donde:
* $I_\alpha$ es la corriente de entrada a la población $\alpha$
* ```math$\mathbf{1}_{N_\alpha}$``` es el vector $N_\alpha \times 1$ de todos unos
* $\alpha \in \{E, I\}$ (excitatoria o inhibitoria)

<image src="branchs.png" alt="Modelo de circuitos neuronales">

<image src="newplot.png" alt="Modelo de circuitos neuronales">