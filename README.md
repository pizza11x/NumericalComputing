# Numerical Calculation University HomeWork 
Track:
1. - Define a row vector T of 30 positive, pseudo-random integers, of 3 digits. (Preliminarily assign a pre-set seed of 6 digits, of which at least 3 distinct, to the function rng. Example: >> rng (143362); do not use 143362 and don't use consecutive digit sequences like 123456, 456789, 765432, etc)

2. - Name X the vector obtained by ordering the vector T and eliminating it then any repetitions and all components multiples of 6.
- Define m and M respectively as the minimum and maximum of X.
- (check that X has between 20 and 30 components. Otherwise assign another value at the seed in point 1.)

3. - Define a handle to the function f = g (p), composed of a function trigonometric g and a polynomial p.
To define f fix:
• p as a second degree polynomial with coefficients other than 0, ± 1;
• g as one of the functions: sine, cosine, tangent and cotangent.
Example: g (x) = cot (x), p (x) = 3x2 - 7x - 4 and f (x) = cot (3x2 - 7x - 4).
N.B. use a function other than the example.
N.B. it is not necessary to explicitly define p and g
(Verify that f changes sign at least once in X. Evaluate for this purpose
the signs of the components of the vector Y which contains the values corresponding to X
through f).

4. - Find a zero of f in [m, M] with a suitable numerical method (to describe). Set the accuracy as desired between 10−4 and 10−8.

5. - Graph the function f in [m, M].
- Highlight the points of f of the abscissa in the vector Y at point 2.
- Also highlight the zero of the function.

6. - Set at will two non-integer values L and R, with L <3.4 and R> 10.5.
- Generate a square matrix A of order 5, pseudo-random with entries uniformly distributed in [L, R].
- Check that A is not singular, otherwise change the again value of the seed in rng immediately before generating A.

7. - Generate a column vector, to be named x true, of 5 components whole in an interval of your choice.
- Determine the column vector b so that x true is the solution of linear system Ax true = b

8. - Call x sol the solution of the system Ax = b calculated with one or more of the methods studied (Briefly describe the method (s) in the comments)
- Calculate the error that is committed by taking x sol in as the solution place of x true.

9. - Invert two lines of A.
- Replace the income of A below average with the average value (average denotes the mean between L and R).

10. - Define x as a vector extracted from X by choosing ten components non-consecutive of the latter (set the 10 components as desired, or pseudo-casually).
- Define the vector xx as a dense grid in [min (x), max (x)] of N points (fix N> 120).

11. - Interpolate f in x with a polynomial pol int of suitable degree.
- Approximate f in x with a polynomial pol app of degree 3.
- Approximate f in x with a trig polynomial of degree 2.

12. - Draw graphically in a second figure:
• the points of the graph of f in x;
• the points of the graph of pol int in xx;
• the points of the graph of pol app in xx;
• the points of the trig graph in xx.

13. - Define x averaged as the output of the function:
mx = NUMMAT weighted average (x) where the input is the vector x.

14. - Define sxy as the function output:
s = saxpy NUMMAT (x, y, a) where x and y are the two vectors previously defined and a is a number real chosen at random.

15. - Develop the function mx = NUMMAT weighted average (x). (NUMMAT indicates the last 6 digits of your serial number.).
This function receives a vector of n components x = (x1, x2,..., Xn) e outputs an mx vector of components:
mx(i) = (x(1)+x(2))/2             i = 1
        (x(i−1)+2x(i)+x(i+1))/4   2 ≤ i ≤ n − 1
        (x(n−1)+x(n))/2            i = n
        
16. – Sviluppare la funzione s=saxpy NUMMAT(x,y,a). (NUMMAT indica le ultime 6 cifre della propria matricola.).
Tale funzione riceve due vettori di n componenti x = (x1, x2, . . . , xn) e y = (y1, y2, . . . , yn), un numero reale a e restituisce in output un vettore s di componenti:
s(i) = a · x(i) + y(i), i = 1, 2, . . . , n
