# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 23:27:17 2024

@author: Ads
"""

import numpy as np
import matplotlib.pyplot as plt

#Question 1 
def generate_transition_matrix(p,q, taille):
    #taille= nombre d'états / cardinal de l'espace d'état E
    m=taille
    P=np.diag(q*(1-p)*np.ones(m),-1)+np.diag((q*p+(1-p)*(1-q))*np.ones(m+1))+np.diag(p*(1-q)*np.ones(m),1)
    P[0,0]=1-p
    P[0,1]=p
    P[-1,-1]=((q*(1-p))/((q*(1-p)+ q*p+(1-p)*(1-q)))) # pour avoir une matrice de transition, il faut normaliser
    P[-1,m-1]=(q*p+(1-p)*(1-q))/((q*(1-p)+ q*p+(1-p)*(1-q)))

                               
    return P


def simulate_markov_chain(transition_matrix,etat_initial, n):
    
    taille = len(transition_matrix)
    etat = [etat_initial]
    etat_actuel = etat_initial
    
    for i in range(n):
        # Générer la prochaine étape selon la distribution de probabilité spécifiée par la matrice de transition
        prochain_etat = np.random.choice(range(taille), p=transition_matrix[etat_actuel])
        etat.append(prochain_etat)
        etat_actuel = prochain_etat
    
    return etat

# simulation d'une trajectoire de longueur 100 :
    
etat_initial=0
p=1/2
q=3/4
taille=5
n=100

transition_matrix= generate_transition_matrix(p,q,taille)

trajectoire =simulate_markov_chain(transition_matrix, etat_initial, n) 


plt.figure(figsize=(10, 6))
plt.plot(trajectoire,marker='o')
plt.title('Trajectoire de la chaîne de Markov')
plt.xlabel('Étapes')
plt.ylabel('États')
plt.grid(True)
plt.show()

#affichage de la trajectoire sous forme d'une liste

print("la trajectoire de la chaine de markov est :", trajectoire) 

#Méthode 2 pour simuler une trajectoire :
"""def traj(p,q,m,X0):
    #X0 1ère état 
    X=[X0]
    #P matrice de transition
    P=np.diag(q*(1-p)*np.ones(m),-1)+np.diag((q*p+(1-p)*(1-q))*np.ones(m+1))+np.diag(p*(1-q)*np.ones(m),1)
    P[0,0]=1-p
    P[0,1]=p
    for i in range (m):
        X.append(np.random.choice(np.arange(m+1),p=P[X[i],:])) 
    return X """

    
#Question 2

p=1/2
q=3/5

taille = 3
P=  generate_transition_matrix(p,q,taille)
print(P)

B= np.linalg.matrix_power(P,1000000) # si (Xn) est irreductible et apériodique sur E fini
                                     # alors cette puissance tend vers la mesure stationnaire

print(B)

# valeur exacte de la msure invariante 

def mesure_invariante(p,q,k):
    rho= p*(1-q)/(q*(1-p))
    

    pi0=((1-q)*(1-rho)+rho)/((1-q)*(1-rho))
    
    pi1= 1/ pi0
    
    if k==0:
        return pi1
    else :
        return (1/(1-q))*(rho**k)*pi1
    
mesure_invariante(p, q, 2)

print (mesure_invariante(p, q, 2))



"""
code de deucresefand 
def proba_stat(p,q,k):
    pbar=1-p
    qbar=1-q
    rho=p*qbar/pbar/q
    if k==0:
        return 1-p/q
    else:
        return (1-p/q)/qbar*rho**k

proba_stat(1/2,3/5,2)  """

#Question 3


k=3
p=1/2
q=3/5


# Linalg ne marche pas!
A = P - np.ones(4)

b = np.zeros(4)

x = np.linalg.solve(A,b)

print (x)

#calcul de mesure invariante exacte

def proba_stat_tronquee(p,q):
    pbar=1-p
    qbar=1-q
    rho=p*qbar/pbar/q
    cumul=1
    for j in range(1,4):
        cumul+=rho**j/qbar
    a=[1/cumul]   
    for j in range(1,4):
        a.append(1/qbar/cumul*rho**j)
    return a

C= proba_stat_tronquee(1/2,3/5)  
print(C)

# question 3.3)

n= 1000 

def trajectoire(etat, n):
    return simulate_markov_chain(P,etat, n)


L=[trajectoire(i, n).count(i)/n for i in range (4)] 
print(L)

#question 3.4

for n in [100, 500, 1000, 10000, 100000]:
    
     
     L1=[trajectoire(i, n).count(i)/n for i in range (4)] 
     print(L1)
     L2= C
     x=[1,2,3,4]
     print(L2)    
     plt.plot(x, L1,label=f'valeur approchée pour n={n}'  )
     plt.plot(x,L2)
     plt.legend() 
     
plt.show()


print(C)

#Temps moyen de séjour dans un état donné

for i in range (4):
    
    print ('le temps moyen de sejour dans l"état ',i, 'est :', 1/(1-P[i,i]))
     
    
    
    
    
    
    
    
    
    
