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

# - basis()

## l'implementation de lalgorithme de "cox boor" pour generer les valeurs des fonctions de base
 c        = order of the B-spline basis function
 
 d        = first term of the basis function recursion relation
 
 e        = second term of the basis function recursion relation
 
 npts     = number of defining polygon vertices
 
 n[]      = array containing the basis functions
 
 n[1] contains the basis function associated with B1 etc.
       
 nplusc   = constant -- npts + c -- maximum number of knot values
 
 t        = parameter value
 
 temp[]   = temporary array
 
 x[]      = knot vector
 
 ```c
 void basis(int c, float t, int npts, int x[], float n[])
{
    int nplusc;
    int i, k;
    float d, e;
    float temp[20000];
    nplusc = npts + c;
    /* calculate the first order basis functions n[i][1]    */
    for (i = 1; i <= nplusc - 1; i++)
    {
        if ((t >= x[i]) && (t < x[i + 1]))
            temp[i] = 1;
        else
            temp[i] = 0;
    }
    /* calculate the higher order basis functions */
    for (k = 2; k <= c; k++)
    {
        for (i = 1; i <= nplusc - k; i++)
        {
            if (temp[i] != 0) /* if the lower order basis function is zero skip the                         calculation */
                d = ((t - x[i]) * temp[i]) / (x[i + k - 1] - x[i]);
            else
                d = 0;
            if (temp[i + 1] != 0) /* if the lower order basis function is zero skip the calculation  */
                e = ((x[i + k] - t) * temp[i + 1]) / (x[i + k] - x[i + 1]);
            else
                e = 0;

            temp[i] = d + e;
        }
    }

    if (t == (float)x[nplusc])
    { /*    pick up last point    */
        temp[npts] = 1;
    }

    /* put in n array   */

    for (i = 1; i <= npts; i++)
    {
        n[i] = temp[i];
    }
}
 ```
 
 

# - bspline

## 


# Comparaison





# performance
