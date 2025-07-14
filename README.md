# Universidad Nacional de Colombia

## Facultad de Ciencias

## Programa de Ciencias de la Computacion

# Análisis de Bifurcaciones

## Presentado por: **Jacobo Ayala Giraldo**

## Materia: **Modelado y Simulación**

---

# 1. Introducción y Problema Central

El cerebro se organiza a distintas escalas, siendo la **mesoscópica** (circuitos de decenas a miles de neuronas) fundamental para conectar la neurona individual con el comportamiento complejo.

Las herramientas matemáticas estándar, como las aproximaciones de **campo medio** (mean-field), son efectivas para redes muy grandes pero fallan al describir circuitos pequeños. Estas aproximaciones pueden ocultar fenómenos importantes que solo ocurren en redes de tamaño finito, los llamados **efectos de tamaño finito** (finite-size effects).

**Objetivo del trabajo:** Desarrollar y utilizar un método matemático para analizar la dinámica de redes neuronales pequeñas, demostrando que su comportamiento es mucho más rico y complejo de lo que predice la teoría de campo medio.

---

# 2. El Modelo Utilizado

Se emplea un **modelo de tasa de disparo** (*firing-rate model*), que describe la actividad promedio de las neuronas.

El circuito se compone de dos poblaciones neuronales: una **Excitatoria (E)** y una **Inhibitoria (I)**.

**Supuestos Clave:**
- **Homogeneidad:** Todas las neuronas dentro de una misma población son idénticas, lo que dota al modelo de una alta simetría.
- **Conectividad Total:** Todas las neuronas están conectadas entre sí.
- **Determinista:** El modelo no considera ruido o fluctuaciones aleatorias.

La siguiente figura muestra el esquema del modelo, con las interacciones entre y dentro de las poblaciones E e I:

<image src="image.png" alt="Modelo de circuitos neuronales">

## Ecuaciones del modelo

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

# 3. El Fenómeno Central: Ruptura Espontánea de Simetría

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

- **Régimen de Inhibición Fuerte:** Se produce una **ruptura espontánea de simetría**: aunque todas las neuronas inhibitorias son idénticas, pueden adoptar diferentes niveles de actividad.

La siguiente figura ofrece una analogía intuitiva: con inhibición fuerte, el "paisaje de energía" del sistema desarrolla múltiples valles, y el sistema "cae" en uno de ellos, rompiendo la simetría original.

<image src="branchs.png" alt="Ramas de soluciones y ruptura de simetría">

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

Como resultado, en el diagrama de bifurcación aparecen **ramas secundarias** de soluciones que emanan de **puntos de ramificación** (BP), representando los nuevos estados asimétricos.

---

# 4. Resultados Funcionales y Relevancia Biológica

## Histéresis y Memoria de Trabajo

La estructura plegada de las soluciones conduce a la **histéresis**, donde el estado de la red depende de su historia reciente. Al aumentar un estímulo ($I_E$), la red salta a un estado de alta actividad, y al disminuirlo, se mantiene en ese estado antes de caer, demostrando un comportamiento con "memoria". Este mecanismo se ha propuesto como una base para la **memoria de trabajo**.

<image src="figura_6_recreada.png" alt="Histéresis y memoria de trabajo">

## Flexibilidad Computacional

La red no tiene un único modo de operación; su comportamiento puede cambiar dinámicamente. La superficie de soluciones, o **colector de catástrofe**, muestra que la red puede actuar como:
- Integrador con fugas (Leaky Integrator)
- Interruptor (Switch)
- Integrador perfecto (Perfect Integrator)

Esta flexibilidad computacional es clave, y se ha sugerido que este tipo de mecanismo podría explicar la **percepción de intervalos de tiempo en el cerebro**.


---

# 5. Conclusiones Finales

Los circuitos neuronales pequeños exhiben **dinámicas complejas y efectos de tamaño finito** que son ignorados por las aproximaciones de campo medio.

La **ruptura espontánea de simetría**, impulsada por una inhibición fuerte, es un fenómeno central que genera heterogeneidad funcional a partir de una estructura homogénea.

Incluso en regímenes más simples, la red muestra comportamientos computacionalmente ricos y biológicamente relevantes como la **histéresis (memoria)** y la **flexibilidad computacional (percepción del tiempo)**.
