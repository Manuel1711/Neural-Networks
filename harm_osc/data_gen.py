#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Definisci la funzione per l'integrando
def integrand(omega, tau, T, rho):
    return (np.exp(-omega * tau) + np.exp(-(T - tau) * omega)) * rho(omega) 

# Definisci la funzione per C(tau)
def correlation(tau, T, rho):
    result, _ = quad(lambda omega: integrand(omega, tau, T, rho), 0, np.inf)
    return result

# Definisci la funzione rho(omega) con parametri personalizzabili
def rho(omega, mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3):
    gaussian1 = a1 * np.exp(-((omega - mu1) / sigma1)**2) / (sigma1 * np.sqrt(2 * np.pi))
    gaussian2 = a2 * np.exp(-((omega - mu2) / sigma2)**2) / (sigma2 * np.sqrt(2 * np.pi))
    gaussian3 = a3 * np.exp(-((omega - mu3) / sigma3)**2) / (sigma3 * np.sqrt(2 * np.pi))
    return gaussian1

# Numero di funzioni rho casuali da generare
num_random_rho = 100#100000

# Valori di tau di interesse da 0 a 199
tau_values = range(200)

# Numero di sottoplot per riga
num_subplot_per_row = 8

# Calcola il numero di righe necessarie per tutti i tau
num_rows = len(tau_values) // num_subplot_per_row + 1

# Array bidimensionale per memorizzare i risultati degli integrali dei correlatori
correlation_integrals = np.zeros((len(tau_values), num_random_rho))

# Matrice per memorizzare i parametri utilizzati
parametri_rho = np.zeros((num_random_rho, 9))  # 9 parametri: mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3

for i in range(num_random_rho):
    #print(i)
    # Genera valori casuali per i parametri della funzione rho
    mu1 = i/num_random_rho
    mu2 = np.random.uniform(0, 1)
    mu3 = np.random.uniform(0, 1)
    #sigma1 = np.random.uniform(0, 0.1)
    #sigma2 = np.random.uniform(0, 0.1)
    #sigma3 = np.random.uniform(0, 0.1)
    sigma1=0.01
    sigma2=0.01
    sigma3=0.01
    a1 = np.random.uniform(0.1, 1.0)
    a2 = np.random.uniform(0.1, 1.0)
    a3 = np.random.uniform(0.1, 1.0)

    # Memorizza i parametri nella matrice parametri_rho
    parametri_rho[i] = [mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3]
    
    for j, tau in enumerate(tau_values):
        # Calcola l'integrale del correlatore per la funzione rho corrente e il valore di tau corrente
        C_integral = correlation(tau, 200, lambda omega: rho(omega, mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3))
        
        # Assegna il risultato all'array bidimensionale
        correlation_integrals[j, i] = C_integral



# Grafica gli istogrammi degli integrali dei correlatori per diversi valori di tau
plt.figure(figsize=(16, 12))
for i in range(len(tau_values)):
    plt.subplot(num_rows, num_subplot_per_row, i + 1)
    plt.hist(correlation_integrals[i], bins=20, edgecolor='k')
    plt.xlabel('Integrale del Correlatore', fontsize=8)
    plt.ylabel('Frequenza', fontsize=8)
    plt.title(f'Tau = {tau_values[i]}', fontsize=10)
    plt.grid(True)

plt.tight_layout()
plt.show()


# In[23]:

'''
# Grafica C(\tau) in funzione di \tau
plt.figure(figsize=(8, 6))
plt.plot(tau_values, correlation_integrals[:,99], label=r'$C(\tau)$', linewidth=2)
plt.xlabel(r'$\tau$', fontsize=14)
plt.ylabel(r'$C(\tau)$', fontsize=14)
plt.title('Funzione di Correlazione $C(\\tau)$', fontsize=16)
plt.grid(True)
plt.legend()
plt.show()
'''

# In[24]:


for i in range(num_random_rho):

    correlation_values = correlation_integrals[:, i]
    omega_values = np.linspace(0, 1, 1000)
    mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3 = parametri_rho[i]

    # Apri il file per la scrittura
    fileop = f'fdata/Data_{i}_1_reg.txt'
    with open(fileop, 'w') as file:
        # Scrivi i primi 100 valori di correlazione
        for value in correlation_values:
            file.write(f"{value}\n")
    
        # Aggiungi una riga vuota
        file.write("\n")
    
        # Scrivi i valori di rho(omega) per i primi 100 valori di omega
        for omega in omega_values:
            rho_value = rho(omega, mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3)
            file.write(f"{rho_value}\n")

'''
mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3 = parametri_rho[10]

omega_values = np.linspace(0, 1, 1000) 
rho_values=rho(omega_values, mu1, mu2, mu3, sigma1, sigma2, sigma3, a1, a2, a3)
# Grafica \rho(\omega) in funzione di omega
plt.figure(figsize=(8, 6))
plt.plot(omega_values, rho_values, label=r'$\rho(\omega)$', linewidth=2)
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\rho(\omega)$', fontsize=14)
plt.title('Densità Spettrale $\\rho(\\omega)$', fontsize=16)
plt.grid(True)
plt.legend()
#plt.show()
'''

# In[ ]:




