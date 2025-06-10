import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

class MesoscopicNeuralNetwork:
    def __init__(self, J_II=-10, I_E= 10, I_I = -10):
        # Parámetros del modelo (basados en el paper)
        self.N_E = 8  # Poblaciones excitatorias
        self.N_I = 2  # Poblaciones inhibitorias
        self.N = self.N_E + self.N_I

        # Constantes de tiempo
        self.tau_E = 1
        self.tau_I = 1

        # Pesos sinápticos
        self.J_EE = 10
        self.J_EI = -70  # Inhibitorio
        self.J_IE = 70
        self.J_II = J_II # Inhibitorio, parametro de bifurcacion

        # Parámetros de activación
        self.nu_max_E = 1
        self.nu_max_I = 1
        self.Lambda_E = 2
        self.Lambda_I = 2
        self.V_t_E = 2
        self.V_t_I = 2

        # Entrada externa
        self.I_E = I_E  #Parametro de bifurcation
        self.I_I = I_I #Parametro de bifurcation

    def activation_function(self, V, nu_max, Lambda, V_t):
        """Función de activación algebraica suave"""
        x = V - V_t
        return (nu_max/2) * (1 + (Lambda/2) * x / np.sqrt(1 + (Lambda**2/4) * x**2))

    def F(self, mu_E, mu_I):
        """Nullcline para población excitatoria"""
        A_E = self.activation_function(mu_E, self.nu_max_E, self.Lambda_E, self.V_t_E)
        A_I = self.activation_function(mu_I, self.nu_max_I, self.Lambda_I, self.V_t_I)

        return (-mu_E/self.tau_E +
                (self.N_E - 1)/(self.N - 1) * self.J_EE * A_E +
                self.N_I/(self.N - 1) * self.J_EI * A_I +
                self.I_E)

    def G(self, mu_E, mu_I):
        """Nullcline para población inhibitoria"""
        A_E = self.activation_function(mu_E, self.nu_max_E, self.Lambda_E, self.V_t_E)
        A_I = self.activation_function(mu_I, self.nu_max_I, self.Lambda_I, self.V_t_I)

        return (-mu_I/self.tau_I +
                self.N_E/(self.N - 1) * self.J_IE * A_E +
                (self.N_I - 1)/(self.N - 1) * self.J_II * A_I +
                self.I_I)

    def find_fixed_points(self):
        """Encuentra puntos fijos del sistema"""
        def system(vars):
            mu_E, mu_I = vars
            return [self.F(mu_E, mu_I), self.G(mu_E, mu_I)]

        # Múltiples condiciones iniciales para encontrar diferentes puntos fijos
        initial_guesses = [
            [0, 0], [5.7, 40], [-10, -10], [20, -5], [-5.8, 40]
        ]

        fixed_points = []
        for guess in initial_guesses:
            try:
                sol = fsolve(system, guess, xtol=1e-12)
                # Verificar que es realmente una solución
                if np.allclose(system(sol), [0, 0], atol=1e-6):
                    # Evitar duplicados
                    is_duplicate = False
                    for fp in fixed_points:
                        if np.allclose(sol, fp, atol=1e-3):
                            is_duplicate = True
                            break
                    if not is_duplicate:
                        fixed_points.append(sol)
            except:
                continue

        return np.array(fixed_points)

    def plot_nullclines_and_fixed_points(self, I_E_value, ax=None):
        """Grafica nullclines y puntos fijos para un valor dado de I_E"""
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))

        # Actualizar parámetro
        self.I_E = I_E_value

        # Rango para las gráficas
        muE_range = np.linspace(-8, 10, 10000)
        muI_range = np.linspace(-15, 60, 10000)

        # Calcular nullclines
        F_values = []
        G_values = []

        for mu_E in muE_range:
            # Para nullcline E: encontrar mu_I tal que F_E(mu_E, mu_I) = 0
            def eq_E(mu_I):
                return self.F(mu_E, mu_I[0])

            try:
                mu_I_sol = fsolve(eq_E, [0], xtol=1e-8, maxfev=1000,)[0]
                # Verificar que la solución es válida
                if abs(self.F(mu_E, mu_I_sol)) < 1e-6:
                    F_values.append([mu_E, mu_I_sol])
            except:
                pass

        for mu_I in muI_range:
            # Para nullcline I: encontrar mu_E tal que F_I(mu_E, mu_I) = 0
            def eq_I(mu_E):
                return self.G(mu_E[0], mu_I)

            try:
                mu_E_sol = fsolve(eq_I, [0], xtol=1e-10)[0]
                if abs(self.G(mu_E_sol, mu_I)) < 1e-6:
                    G_values.append([mu_E_sol, mu_I])
            except:
                pass

        # Convertir a arrays
        if F_values:
            F = np.array(F_values)
            ax.plot(F[:, 0], F[:, 1], 'm-', linewidth=2,
                   label=r'$\mathcal{F}(\mu_E,\mu_I)=0$')

        if G_values:
            G = np.array(G_values)
            ax.plot(G[:, 0], G[:, 1], 'g-', linewidth=2,
                   label=r'$\mathcal{G}(\mu_E,\mu_I)=0$')

        # Encontrar y graficar puntos fijos
        fixed_points = self.find_fixed_points()
        if len(fixed_points) > 0:
            ax.plot(fixed_points[:, 0], fixed_points[:, 1], 'ko',
                   markersize=8, markerfacecolor='black', markeredgecolor='white',
                   markeredgewidth=2, label='Puntos fijos')

        # Configuración de la gráfica
        ax.set_xlabel(r'$\mu_E$', fontsize=14)
        ax.set_ylabel(r'$\mu_I$', fontsize=14)
        ax.set_title(f'$I_E = {I_E_value}$', fontsize=16)
        ax.legend(fontsize=12)
        ax.set_xlim(-8, 10)
        ax.set_ylim(-15, 60)

        return fixed_points

def create_bifurcation_diagram():
    """Crea el diagrama de bifurcación completo"""
    network = MesoscopicNeuralNetwork()
    network_1 = MesoscopicNeuralNetwork(I_E=13)

    # Crear subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax1, ax2 = axes[0]
    ax3, ax4 = axes[1]


    # Primer caso: I_E = 10
    print("Analizando I_E = 10...")
    fp1 = network.plot_nullclines_and_fixed_points(10, ax1)
    print(f"Puntos fijos encontrados: {len(fp1)}")
    if len(fp1) > 0:
        for i, fp in enumerate(fp1):
            print(f"  Punto {i+1}: μ_E = {fp[0]:.3f}, μ_I = {fp[1]:.3f}")

    # Segundo caso: I_E = 13
    print("\nAnalizando I_E = 13...")
    fp2 = network_1.plot_nullclines_and_fixed_points(13, ax2)
    print(f"Puntos fijos encontrados: {len(fp2)}")
    if len(fp2) > 0:
        for i, fp in enumerate(fp2):
            print(f"  Punto {i+1}: μ_E = {fp[0]:.3f}, μ_I = {fp[1]:.3f}")

    plt.tight_layout()
    plt.suptitle('Análisis de Bifurcaciones en Redes Neuronales Mesoscópicas',
                 fontsize=18, y=1.02)
    
    IE_vector = np.linspace(-4, 18, 10000)
    all_fixed_points = []

    for i in IE_vector:
        network.I_E = i
        fixed_points = network.find_fixed_points()
        for fp in fixed_points:
            all_fixed_points.append([i,fp[0],fp[1]])

    all_fixed_points = np.array(all_fixed_points)

    
    ax3.scatter(all_fixed_points[:,0], all_fixed_points[:,1], color='black',alpha=0.3, label=r'$\mu_E$')
    ax3.set_ylabel(r'$\mu_E$', fontsize=14)
    ax3.legend(fontsize=12)
    ax3.set_xlim(-4, 18)
    ax3.set_ylim(-4, 10)

    ax4.scatter(all_fixed_points[:,0], all_fixed_points[:,2], color='black',alpha=0.3, label=r'$\mu_I$')
    ax4.set_xlabel(r'$I_E$', fontsize=14)
    ax4.set_ylabel(r'$\mu_I$', fontsize=14)
    ax4.legend(fontsize=12)
    ax4.set_xlim(-4, 18)
    ax4.set_ylim(-15, 60)
    plt.show()

    return fig

# Ejecutar análisis
if __name__ == "__main__":
    # Crear diagrama principal
    fig = create_bifurcation_diagram()