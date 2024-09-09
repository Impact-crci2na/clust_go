from collections import Counter
import requests
import json
from goatools import obo_parser
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from bioservices.uniprot import UniProt
u = UniProt(verbose = False)

"""### ▲ FONCTIONS

#### ▬ Conversion d'identifiant
"""

def convert_ACC_GENENAME(id) :
    convert = u.mapping(fr="UniProtKB_AC-ID", to="Gene_Name", query=id)
    return convert['results'][0]['to']

"""#### ▬ Lecture du fichier"""
def import_data(path):
    list_acc=[]
    with open(path,"r") as f:
        for l in f.readlines():
            list_acc.append(l[:-1])
    print("Number of Protein :", len(list_acc))
    return(list_acc)
"""#### ▬ Association des GO terms aux protéines"""

def create_dico_go_term(list_id, type_of_go_terms,GO_file ="../go-basic.obo" ):

    go = obo_parser.GODag(GO_file)
    ## VALUES of the dictionnary
    dico_go_term = {}
    for i in list_id:
        list_temp = []
        r=requests.get(f"http://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId={i}")
        data=r.json()
        try :
            for j in data["results"]:
                if j["goAspect"] == type_of_go_terms:
                    go_name=go[j["goId"]].name
                    list_temp.append(go_name)
                dico_go_term[i]=list_temp
        except :
            print(j)
    ## KEYS of the dictionnary
    list_go=[]
    for v in dico_go_term.values():
        for v1 in v:
            if v1 not in list_go:
                list_go.append(v1)
    #print((len(list_go)))
    return dico_go_term, list_go

def define_matrix(dico_go_term, list_go) :
    list_prot=list(dico_go_term.keys())
    array=np.zeros((len(list_prot),len(list_go)))
    for i in range(len(list_prot)):
        for j in range(len(list_go)):
            if list_go[j] in dico_go_term[list_prot[i]]:
                array[i][j]=1
            else:
                array[i][j]=0

    df = pd.DataFrame(array, index=list_prot,columns=list_go)
    #df.shape
    #print(df.head())
    return df

"""#### ▬ Clustering via la méthode des KMeans"""

def clustering_kmeans(dico_go_term, list_go, df,nb_cluster= 2,top = 5, file_name = "cluster_proteine_go_terms") :
    print(nb_cluster)
    kmeans = KMeans(n_clusters=nb_cluster, random_state=0).fit(df)
    pca = PCA(n_components=2)
    data=pca.fit_transform(df)

    df_go = df
    df_go["cluster"]=kmeans.labels_
    colormap=np.array(['Red','green','blue','orange','black','purple','brown','gray','olive','cyan'])
    res=kmeans.__dict__
    with open(f"{file_name}.txt","w") as f:
      for k in range(nb_cluster):
        print(k)
        f.write("\n\n")
        f.write("*"*50)
        f.write(f"\nClasse: {k} ({colormap[k]})\n")
        f.write("*"*50)
        f.write("\n\n")
        f.write("\nListe des keywords:")
        list_go_weight=[]
        for i in np.argsort(res['cluster_centers_'][k]*-1):
            list_go_weight.append(list_go[i])
        for w in list_go_weight[:top]:
            f.write(f"\n{w}")
        f.write(f"\n\nListe des protéines de la classe {k}:\n")
        for p in df_go.index[df_go["cluster"]==k].tolist():
          f.write(f"{p} ")
    ## Résultats
    plt.scatter(data[:,0],data[:,1],c=colormap[kmeans.labels_])
    plt.savefig(f"{file_name}_ACP_go_terms.png", format="PNG")
    data_to_write = []

 ###### il faudra rajouter la partie avec l'ériture du fichier json .
    return df_go

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

def clustering_kmeans_best(df):
    # Plage du nombre de clusters à tester
    range_n_clusters = range(2, 11)  # tester de 2 à 10 clusters
    # Initialisation des listes pour stocker les résultats des scores
    silhouette_scores = []
    #sse = []  # Sum of Squared Errors (SSE)
    best_score = -1
    best_n_clusters = 0
    for n_clusters in range_n_clusters:
        print(f"Processing {n_clusters} clusters...")
        if n_clusters >= df.shape[0] :
            print("invalid number of cluster")
        else:
            kmeans = KMeans(n_clusters=n_clusters, random_state=0)
            cluster_labels = kmeans.fit_predict(df)
        # Calcul du score de silhouette pour l'évaluation de la qualité du clustering
            silhouette_avg = silhouette_score(df, cluster_labels)
            silhouette_scores.append(silhouette_avg)
        if silhouette_avg > best_score:
            best_score = silhouette_avg
            best_n_clusters = n_clusters
        # Calcul de la SSE pour la méthode du coude
        #sse.append(kmeans.inertia_)  # 'inertia_' est la SSE
        # calcule du meilleurs score
    print(f"Le meilleur score de silhouette est {best_score} pour {best_n_clusters} clusters.")

    # Affichage du graphique du score de silhouette
    # plt.figure(figsize=(10, 5))
    # plt.plot(range_n_clusters, silhouette_scores, marker='o')
    # plt.title('Scores de silhouette pour différents nombres de clusters')
    # plt.xlabel('Nombre de Clusters')
    # plt.ylabel('Score de Silhouette')
    # plt.grid(True)
    # plt.show()
    return best_n_clusters


def take_top(data):
    # Compter le nombre de membres par cluster
    cluster_counts = data['cluster'].value_counts()
    print("cluster count")
    # Identifier le cluster avec le plus grand nombre de membres
    max_cluster = cluster_counts.idxmax()
    # Récupérer les membres du cluster le plus grand
    largest_cluster_proteins = data.index[data['cluster'] == max_cluster].tolist()
    print("largest")
    print(largest_cluster_proteins)
    print("max cluster")
    print(max_cluster)
    return largest_cluster_proteins

def determine_highlight(data_list, name):
  """
  Analyzes a list and highlights elements that appear more than once.

  Args:
      data_list: The list of data to analyze.
      name: A name to use for file creation (optional).

  Returns:
      None
  """

  try:
    # Print the original list (optional)
    # print(data_list)
    data_list_second=[]
    for i in data_list:
        data_list_second+=i
    # Create a frequency count of elements
    element_counts = Counter(data_list_second)

    # Find elements that appear more than once
    highlighted_elements = [item for item, count in element_counts.items() if count > 1]

    # Handle file creation (optional):
    if name:
      create_list_of_list(highlighted_elements, f"{name}")
      print("creation of files")

  except Exception as e:
    print(f"Error: {e}")





def create_list_of_list(data, file_name = "top"):
    try:
        with open(f"{file_name}.txt","w") as f:
            for prot in data:
                f.write(prot)
                f.write("\n")
    except:
        print(f"an error had occured")
