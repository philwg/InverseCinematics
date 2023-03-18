#================================================
# Un fichier Makefile permet d'automatiser la compilation
# d'un projet. Il fonctionne avec un systeme de Regle
# et Dependance
#
# Regle: Dependance
# 	 Commande
#
# Donc attention la Dependance est appellé avant 
# d'effectuer la commande
#
# Chaque regle peut être appellée independament 
# Un makefile peut être fragmenté, avoir des conditions
# ect... cependant pour les petits projet ce makefile suffit
#
# Un bon tuto ce trouve sur developpez.com
# "Introduction à Makefile" de gl
#
# Pour executer un makefile:
# make <Nom de la Regle>
# Si aucune Regle est precisée la regle all est
# celle par defaut
#
# Si apres compilation tu ne change qu'un header
# appelle la regle rebuild avec make rebuild
# en effet normalement on peut faire dependre
# les sources de leurs headers pour que make build
# a nouveau le projet si un header est changé mais
# j'ai eu la flemme... Bonne chance est hesite pas
# à demandé de l'aide
#================================================

#================================================
# Commande du compilateur
# gcc => C
# g++ => C++
# javac => Java
# ...
#================================================
CC=g++

#================================================
# Flags de compilation
# -Wall => Warning All (Affiche tout les warning)
#================================================
CFLAGS=-Wall -g

#================================================
# Flags de link
# -lGL => OpenGL
# -lGLU => GLU
# -lglut => freeGLUT
#================================================
LDFLAGS=-lGL -lGLU -lglut -O2 -larmadillo -llapack -lblas

#================================================
# Nom du fichier executable final
# Chemin du dossier des sources
#================================================
EXEC=Launch
SRC_DIR=.

#================================================
# Liste des sources (.c)
# Liste des objets (.o)
#================================================
SRC=$(wildcard $(SRC_DIR)/*.cpp)
OBJ=$(SRC:.cpp=.o)

#================================================
# Regle par defaut (make)
# all => Regle all
# $(EXEC) => Appelle de la dependance $(EXEC)
#================================================
all: $(EXEC)

#================================================
# Regle de construction de l'executable
# $(EXEC) => Cible de la regle
# $(OBJ) => Appelle de la dependance $(OBJ) (%.o)
#
# @ Rend makefile silencieux mais pas le compilateur
# $@ Le nom de la cible, fichier à créer (ici $(EXEC) )
# $^ La liste des dépendance (resultat de $(OBJ))
#================================================
$(EXEC): $(OBJ)
	@$(CC) -o $@ $^ $(LDFLAGS) 

#================================================
# $< Le nom de la premiere dependance pour chaque cible
#================================================
%.o: %.c
	@$(CC) -o $@ -c $< $(CFLAGS)

#================================================
# .PHONY Reconstruit systematiquement ses dependances
#================================================
.PHONY: clean mrproper rebuild

#================================================
# Supprime tout les .o
#================================================
clean:
		@rm -rf $(SRC_DIR)/*.o

#================================================
# Appelle la dependance clean et supprime l'executable
#================================================
mrproper: clean
		@rm -rf $(EXEC)

#================================================
# Appelle la dependance clean et all
#================================================
rebuild: clean all