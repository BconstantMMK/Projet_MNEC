# Compilateur utilisé
CC=g++

run : main.cc fonctions1D.cpp eulerpoisson.cpp fonctions1D.h eulerpoisson.h
	$(CC) -O3 -std=c++11  main.cc fonctions1D.cpp eulerpoisson.cpp  -o run

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ run
