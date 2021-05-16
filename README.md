# HiddenMarkovModel

CODE : 

———————————————————————————————————— 

1. backward_forward_discret : 

Fonctions principales : 

backward_prob : renvoie la matrice backward 
forward_prob : renvoie la matrice forward
prob : renvoie la probabilité d’observation d’une séquence à partir des matrices backward et forward et d’un temps t. 

Output du fichier : 

- matrice alpha_t 
- matrice beta_t 
- probabilité d’exécution de la séquence 

———————————————————————————————————— 
2. backward_forward_continu : 
Utilisation du package HiddenMarkov

Fonctions principales : 

backwardforward : renvoie les matrices alpha et beta 

Output du fichier : 

- matrice alpha_t 
- matrice beta_t 
- probabilité d’exécution de la séquence  au temps t = 3

———————————————————————————————————— 

3. EM_discret : 

Fonctions principales : 

backward_prob : renvoie la matrice backward 
forward_prob : renvoie la matrice forward
prob : renvoie la probabilité d’observation d’une séquence à partir des matrices backward et forward et d’un temps t. 
Baum_Welch : estime les paramètres selon les formules de Baum-Welch. 


Output du fichier : 

- Evolution de la proba d’observation : Figure 2 
 

————————————————————————————————————

4. EM_continu : 

Fonctions principales : 

- proba : renvoie une liste contenant les matrice backward et forward

- bwcontrol : fonction qui permet de paramétrer la fonction BaumWelch

- BaumWelch : fonction du package HiddenMarkov qui permet d’exécuter EM sur un modèle HMM

- parameter_list : renvoie l’évolution des paramètres considérés (ici les moyennes des gaussiennes) en fonction des itérations. 
Attention : ici le vecteur de moyenne est (-2,0,2) (conformément à l’article étudié) et n’est pas rentré en paramètre. De même pour la variance, de 0.5. 

- likelihood_list : évolution de la log-vraissemblance complète au cours des itérations. 
Attention : ici le vecteur de moyenne est (-2,0,2) (conformément à l’article étudié) et n’est pas rentré en paramètre. De même pour la variance, de 0.5. 

bootstrap_IC : effectue K estimations des paramètres via EM. 


Output du fichier : 

- Evolutions des paramètre (m_1,m_2,m_3)  (Figure4)
- Evolution de la log-vraissemblance au cours des itérations (Figure 3)
- Histogramme de l’échantillon Bootstrap de m_3 et a_22 Figure 5 (Figure 5 et Figure 6) 
- Intervalle de confiance m_3 et a_22  

—————————————————————————————————————

5. Gibbs :
Attention pour une séquence de longueur t = 1000 et 2000 itérations, le calcul peut être long (environ 1h30). 

Fonctions principales : 

- gene_markov_chain : simulation d’une chaîne de Markov à matrice de transition non-homogène. 
- Gibbs_sampler : échantillonnage de Gibbs selon ordre considérés dans l’article de Tobias Ryden (2008). 

Output du fichier : 

- Graphique de l’évolution du paramètre m_2 (Figure 7). 

—————————————————————————————————————

6. Viterbi :

Fonctions principales : 

- viterbi : calcul la séquence cachée optimale étant donnée une séquence observée ainsi que la probabilité d’observation de cette séquence observée. 

Output du fichier : 

- chemin optimal et probabilité d’observation
