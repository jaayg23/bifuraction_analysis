# Universidad Nacional de Colombia

## Facultad de Ciencias

## Programa de Ciencias de la Computacion

# Análisis de Bifurcaciones en Circuitos Neuronales Mesoscopicos

## Presentado por: **Jacobo Ayala Giraldo**

## Materia: **Modelado y Simulación**

---

## 0. Resumen

Este trabajo investiga la dinámica de redes neuronales de tamaño finito, un área donde las aproximaciones de campo medio tradicionales resultan insuficientes. Se emplea un modelo de tasa de disparo determinista, compuesto por una población excitatoria (E) y una inhibitoria (I), para estudiar cómo la interacción entre neuronas genera comportamientos complejos. Mediante el análisis de bifurcaciones, se demuestra que una fuerte auto-inhibición dentro de la población I provoca una ruptura espontánea de la simetría, dando lugar a nuevos estados de equilibrio asimétricos que no son predichos por la teoría de campo medio. Estos estados generan fenómenos funcionalmente relevantes, como la histéresis, que puede servir como base para la memoria de trabajo, y una notable flexibilidad computacional, permitiendo que el circuito actúe como integrador, interruptor o memoria. Los resultados subrayan la importancia de los efectos de tamaño finito en la comprensión de la computación neuronal.

**Palabras Clave:** Redes Neuronales, Análisis de Bifurcaciones, Ruptura de Simetría, Efectos de Tamaño Finito, Modelo de Tasa de Disparo, Histéresis.

---

## 1. Introducción

El cerebro humano opera a través de una compleja organización jerárquica. La escala **mesoscópica**, que abarca circuitos de decenas a miles de neuronas, es de vital importancia, ya que sirve de puente entre la actividad de neuronas individuales y los comportamientos cognitivos complejos [1]. Históricamente, el análisis de grandes redes neuronales se ha apoyado en aproximaciones de **campo medio** (*mean-field*), las cuales asumen un número infinito de neuronas para simplificar las interacciones [2, 3]. Si bien estas teorías son poderosas, fallan al describir la dinámica de circuitos neuronales pequeños, ocultando fenómenos emergentes que dependen críticamente del número de elementos en la red, conocidos como **efectos de tamaño finito** (*finite-size effects*) [4, 5].

El planteamiento del problema se centra en esta brecha: las herramientas matemáticas estándar son insuficientes para capturar la riqueza dinámica de los circuitos mesoscópicos. El **objetivo principal** de este trabajo es desarrollar y aplicar un método de análisis basado en la teoría de bifurcaciones para demostrar que la dinámica de redes neuronales pequeñas es significativamente más compleja de lo que predice la teoría de campo medio.

El presente artículo se estructura de la siguiente manera: la Sección 2 detalla el modelo matemático y la metodología de análisis. La Sección 3 presenta los resultados centrales, enfocándose en el fenómeno de ruptura de simetría y sus consecuencias funcionales. Finalmente, la Sección 4 expone las conclusiones, limitaciones y posibles líneas de investigación futura.

---

## 2. Planteamiento del Problema y Metodología

### 2.1. Modelo del Circuito Neuronal

Para investigar el fenómeno, se utiliza un **modelo de tasa de disparo** (*firing-rate model*), un enfoque determinista que describe la actividad promedio de poblaciones neuronales [6]. El circuito base consta de dos poblaciones: una **Excitatoria (E)** y una **Inhibitoria (I)**, bajo los siguientes supuestos:

- **Homogeneidad:** Todas las neuronas de una misma población son idénticas.
- **Conectividad Total:** Cada neurona está conectada con todas las demás neuronas de la red.


<image src="image.png" alt="Modelo de circuitos neuronales">

Figura 1. Esquema del modelo de circuito neuronal con poblaciones excitatoria (E) e inhibitoria (I) y sus conexiones sinápticas.

## Ecuaciones del modelo

La dinámica del potencial de membrana $V_i(t)$ de cada neurona $i$ se describe mediante el siguiente sistema de ecuaciones diferenciales:

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

La función de activación utilizada es suave y sigmoidal, definida como:

```math
\phi_{j}(V) = \frac{v_{j}^{max}}{2} \left[ 1 + \frac{\frac{\Lambda_{j}}{2}(V - V_{j}^{T})}{\sqrt{1 + \frac{\Lambda_{j}^{2}}{4}(V - V_{j}^{T})^{2}}} \right]
```

Donde:
* $v_{j}^{max}$ es el valor máximo de la función de activación para la neurona $j$
* $V_{j}^{T}$ es el umbral de activación de la neurona $j$
* $\Lambda_{j}$ es el parámetro de pendiente que controla qué tan abrupta es la transición de la función de activación alrededor del umbral

La matriz de conectividad sináptica se estructura como sigue:

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
* $\mathbf{1}_{N_\alpha}$ es el vector $N_\alpha \times 1$ de todos unos
* $\alpha \in \{E, I\}$ (excitatoria o inhibitoria)

---

## 3. El Fenómeno Central: Ruptura Espontánea de Simetría

El comportamiento de la red depende críticamente de la fuerza de la auto-inhibición ($J_{II}$).

- **Régimen de Inhibición Débil:** Todas las neuronas inhibitorias se comportan de manera idéntica, manteniendo la simetría del sistema. Los estados de equilibrio se encuentran en las intersecciones de las nulclinas, como se muestra en la siguiente figura:

```math
\begin{cases}
\mathcal{F}(\mu_E, \mu_I) \overset{\text{def}}{=} -\frac{1}{\tau_E}\mu_E + \frac{N_E-1}{N-1} J_{EE} \mathcal{A}_E(\mu_E) + \frac{N_I}{N-1} J_{EI} \mathcal{A}_I(\mu_I) + I_E = 0 \\
\\
\mathcal{G}(\mu_E, \mu_I) \overset{\text{def}}{=} -\frac{1}{\tau_I}\mu_I + \frac{N_E}{N-1} J_{IE} \mathcal{A}_E(\mu_E) + \frac{N_I-1}{N-1} J_{II} \mathcal{A}_I(\mu_I) + I_I = 0 .
\end{cases}
```

<image src="fig_2.png" alt="Nulclinas y equilibrio">

Figura 2. Intersección de nulclinas en el espacio de estados, mostrando los puntos de equilibrio simétricos de la red.

- **Régimen de Inhibición Fuerte:** Se produce una **ruptura espontánea de simetría**: aunque todas las neuronas inhibitorias son idénticas, pueden adoptar diferentes niveles de actividad.

La siguiente figura ofrece una analogía intuitiva: con inhibición fuerte, el "paisaje de energía" del sistema desarrolla múltiples valles, y el sistema "cae" en uno de ellos, rompiendo la simetría original.

<image src="branchs.png" alt="Ramas de soluciones y ruptura de simetría">

Figura 3. Diagrama de bifurcación que ilustra la aparición de ramas secundarias y ruptura espontánea de simetría en la red neuronal.

Este fenómeno se visualiza matemáticamente en la siguiente figura. En lugar de simples nulclinas, el sistema se describe por **nulsuperficies** en 3D, cuyas intersecciones revelan nuevos puntos de equilibrio.

```math
\begin{cases}
\mathcal{F}(\mu_E, \mu_{I,0}, \mu_{I,1}) \overset{\text{def}}{=} -\frac{1}{\tau_E}\mu_E + \frac{N_E-1}{N-1} J_{EE} \mathcal{A}_E(\mu_E) + \frac{J_{EI}}{N-1} [ \mathcal{A}_I(\mu_{I,0}) + \mathcal{A}_I(\mu_{I,1}) ] + I_E = 0 \\
\\
\mathcal{G}(\mu_E, \mu_{I,0}, \mu_{I,1}) \overset{\text{def}}{=} -\frac{1}{\tau_I}\mu_{I,0} + \frac{N_E}{N-1} J_{IE} \mathcal{A}_E(\mu_E) + \frac{J_{II}}{N-1} \mathcal{A}_I(\mu_{I,1}) + I_I = 0 \\
\\
\mathcal{H}(\mu_E, \mu_{I,0}, \mu_{I,1}) \overset{\text{def}}{=} -\frac{1}{\tau_I}\mu_{I,1} + \frac{N_E}{N-1} J_{IE} \mathcal{A}_E(\mu_E) + \frac{J_{II}}{N-1} \mathcal{A}_I(\mu_{I,0}) + I_I = 0 .
\end{cases}
\quad (12)
```

<image src="Figura4.png" alt="Nulsuperficies y puntos de equilibrio">

Figura 4. Intersección de nulsuperficies en el espacio de estados tridimensional, revelando los nuevos puntos de equilibrio asimétricos (puntos magenta) que coexisten con el equilibrio simétrico original.

Como resultado, en el diagrama de bifurcación aparecen **ramas secundarias** de soluciones que emanan de **puntos de ramificación** (BP), representando los nuevos estados asimétricos.

---

## 4. Resultados Funcionales y Relevancia Biológica

### Histéresis y Memoria de Trabajo

La estructura plegada de las soluciones conduce a la **histéresis**, donde el estado de la red depende de su historia reciente. Al aumentar un estímulo ($I_E$), la red salta a un estado de alta actividad, y al disminuirlo, se mantiene en ese estado antes de caer, demostrando un comportamiento con "memoria". Este mecanismo se ha propuesto como una base para la **memoria de trabajo**.

<image src="figura_6_recreada.png" alt="Histéresis y memoria de trabajo">

Figura 5. Curva de histéresis que muestra el comportamiento de memoria de trabajo: la red mantiene su estado de alta actividad incluso cuando el estímulo disminuye, hasta un umbral crítico.

### Flexibilidad Computacional

La red no tiene un único modo de operación; su comportamiento puede cambiar dinámicamente. La superficie de soluciones, o **colector de catástrofe**, muestra que la red puede actuar como:
- Integrador con fugas (Leaky Integrator)
- Interruptor (Switch)
- Integrador perfecto (Perfect Integrator)

Esta flexibilidad computacional es clave, y se ha sugerido que este tipo de mecanismo podría explicar la **percepción de intervalos de tiempo en el cerebro**.


---

## 5. Conclusiones Finales

Los circuitos neuronales pequeños exhiben **dinámicas complejas y efectos de tamaño finito** que son ignorados por las aproximaciones de campo medio.

La **ruptura espontánea de simetría**, impulsada por una inhibición fuerte, es un fenómeno central que genera heterogeneidad funcional a partir de una estructura homogénea.

Incluso en regímenes más simples, la red muestra comportamientos computacionalmente ricos y biológicamente relevantes como la **histéresis (memoria)** y la **flexibilidad computacional (percepción del tiempo)**.

## 6. Referencias

[1] R. D. Beer, "On the dynamics of small continuous-time recurrent neural networks," *Adapt. Behav.*, vol. 3, no. 4, pp. 469–509, 1995.

[2] S.-I. Amari, "Dynamics of pattern formation in lateral-inhibition type neural fields," *Biol. Cybern.*, vol. 27, pp. 77–87, 1977.

[3] O. Faugeras, J. Touboul, and B. Cessac, "A constructive mean-field analysis of multi-population neural networks with random synaptic weights and stochastic inputs," *Front. Comput. Neurosci.*, vol. 3, p. 1, 2009.

[4] M. A. Buice and C. C. Chow, "Dynamic finite size effects in spiking neural networks," *PLoS Comput. Biol.*, vol. 9, no. 1, p. e1002872, 2013.

[5] P. C. Bressloff, "Stochastic neural field theory and the system-size expansion," *SIAM J. Appl. Math.*, vol. 70, no. 5, pp. 1488–1521, 2010.

[6] H. R. Wilson and J. D. Cowan, "Excitatory and inhibitory interactions in localized populations of model neurons," *Biophys. J.*, vol. 12, no. 1, pp. 1–24, 1972.

[7] Y. A. Kuznetsov, *Elements of applied bifurcation theory*, vol. 112. New York: Springer-Verlag, 1998.

[8] S. H. Strogatz, *Nonlinear dynamics and chaos*. Sarat Book House, 1994.

[9] G. Deco, et al., "How Local Excitation-inhibition ratio impacts the whole brain dynamics," *J. Neurosci.*, vol. 34, pp. 7886–7898, 2014.

[10] R. D. Beer, "Parameter space structure of continuous-time recurrent neural networks," *Neural Comput.*, vol. 18, pp. 3009–3051, 2006.

[11] D. Hansel and H. Sompolinsky, "Modeling feature selectivity in local cortical circuits," Chap. 13, MIT Press, 1998.

[12] P. C. Bressloff, "Path-integral methods for analyzing the effects of fluctuations in stochastic hybrid neural networks," *J. Math. Neurosci.*, vol. 5, p. 4, 2015.

[13] F. Grimbert, "Mesoscopic models of cortical structures," PhD thesis, Univ. of Nice-Sophia Antipolis, 2008.

[14] O. Faugeras and J. MacLaurin, "Asymptotic description of neural networks with correlated synaptic weights," *Entropy*, vol. 17, pp. 4701–4743, 2015.

[15] B. Cessac, "Increase in complexity in random neural networks," *J. Phys. I (France)*, vol. 5, pp. 409–432, 1995.
