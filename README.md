# B-spline Rationnelles Non uniformes et rendu

# sommaire

- [introduction](#introduction)

- [principe du B-spline](#principe-du-b-spline)
 
- [technique d implementation](#technique-d-implementation)
 
- [Comparaison](#comparaison)
 
- [performance](#performance)



# introduction

La spline de base rationnelle non uniforme ( NURBS ) est un modèle mathématique utilisant des splines de base (B-splines) qui est couramment utilisé en infographie pour représenter des courbes et des surfaces . Il offre une grande flexibilité et précision pour gérer à la fois les formes analytiques (définies par des formules mathématiques courantes ) et modélisées 





# principe du B-spline





# technique d'implementation
- # Knot()
## Knot c'est la fonction qui grnere le knot Vector contenant des valeurs qui seront prises par t selon des intervales donnes

c            = order of the basis function

n            = the number of defining polygon vertices

nplus2       = index of x() for the first occurence of the maximum knot vector value

nplusc       = maximum value of the knot vector -- $n + c$

x()          = array containing the knot vector


```c
void knot(int n, int c, int x[])
{
  int nplusc,nplus2,i;

  nplusc = n + c;
  nplus2 = n + 2;

  x[1] = 0;
  for (i = 2; i <= nplusc; i++){
  if ( (i > c) && (i < nplus2) )
      x[i] = x[i-1] + 1;
  else
     x[i] = x[i-1];
  }
}

```






# Comparaison





# performance
