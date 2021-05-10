clear; clc; close all;
%Antonio Sarno -- 0124001914
%Le function implementate sono:
%-- mx = media_pesata_001914(x)
%-- s = saxpy_001914(x,y,a)

%Punto 1.
%Inizialmente ho aggiornato il valore del seed con la function "rng()",
%In questo caso ho inventato il seed: 357541.
%Dopo ho creato un vettore riga di 30 elementi, i quali sono interi 
%positivi, appartenenti all'intervallo 100 e 999 grazie alla function 
%"randi()".
rng(357541);
T= randi([100 999],1,30);

%Punto 2.
%In questo punto per prima cosa inizializzo un vettore "X" con gli elementi
%del vettore precedentemente creato, utilizzando la function "unique" che
%permette di eliminare gli elementi uguali. Dopo inizia un ciclo for che
%scorre tutta la lunghezza del vettore, e tramite la function "mod()" si
%calcola il resto della divisione tra l'i-esimo elemento e 6 se questo è
%uguale a zero, allora quell'i-esimo elemento viene sostituito con
%un elemento "NaN" (not a number). Finito il ciclo, si eliminano gli elementi
%NaN tramite la function "isfinite()".
%Poi si trovano il minimo e il massimo tramite le funzioni "Min()" e "Max()".
X=unique(T);
for i=1:length(X)
    if(mod(X(1,i),6)==0)
        X(1,i)=NaN;
    end
end
X=X(isfinite(X));
m=min(X);
M=max(X);

%Punto 3.
%Definisco un handle alla funzione "f", che prende in input un solo valore,
%la funzione inventata è "cos((6*x.^2)-3*x-8)".Poi viene assegnato ad un
%vettore riga "Y" i risultati della funzione prendendo in input il vettore "X".
f = @(x) cos((6*x.^2)-3*x-8);
Y = f(X);

%Punto 4.
%Per trovare lo zero di "f" compreso nell'intervallo scelto, utilizzo la
%function di matlab "bisezione", passando come argomenti: la "f", che è la
%nostra funzione; "m" e "M" che sono gli estremi del nostro intervallo; come
%ultimo passiamo l'accuratezza, che ho scelto "1e-6".
%Questo metodo determina il punto medio dell'intervallo, quindi lo divide
%in due sottointervalli di uguale ampiezza, poi sceglie il sottointervallo
%in cui la funzione cambia di segno e continua solo su di esso.
z1= bisezione(f,m,M,1e-6);

%Punto 5.
%Vengono prima di tutto generati 100 punti equispaziati nell'intervallo
%[m,M] tramite la function "linspace()". Successivamente viene trovata la "f" di
%quei punti generati. 
%Di seguito troviamo 3 plot: il primo che disegna graficamente la funzione
%f nei punti generati prima. Il secondo plot evidenzia i punti e il terzo
%evidenzia lo zero della funzione, trovato nel punto4.
x1 = linspace(m,M,100);
y1 = f(x1);
figure
plot(x1,y1,'black');
hold on
grid on
plot (x1, y1,'o');
plot(z1, f(z1), '*-red');
text(z1,0,'Z1');
xlabel('Ascisse');
ylabel('Ordinate');

%Punto 6.
%Definiti i due intervalli "L" e "R", ho generato la matrice quadrata di ordine
%5 con elementi casuali distribuiti uniformemente con la formula 
%" m = a +(b-a).*rand(n) " e dopo ho verificato con la function "det()" 
%il determinante di "A" per vedere se la matrice è non singolare (quindi 
%determinante diverso da zero).
L=3.1; R=11.1;
A=  L + (R-L).*rand(5);
d = det(A);

%Punto 7.
%Generato un vettore colonna, di 5 componenti intere, di un intervallo
%inventato [1,11]. Determiniamo un vettore colonna "b" calcolato con il
%prodotto matriciale righe per colonne tra la matrice quadrata "A" e il
%vettore colonna "x_true".
x_true=randi([1 11],5,1);
b = A*x_true;

%Punto 8.
%Per trovare "x_sol" utilizziamo il metodo di Gauss con pivoting parziale.
%La function "Sgauss()" utilizza a sua volta 2 function al suo interno:
%1) "TriangGauss_pivot()", questa function determina la matrice triangolare
%superiore tramite la tecnica del pivoting; la tecnica del pivoting esamina
%la colonna 'k' della matrice, e identifica l'elemento massimo (in valore
%assoluto) di quella porzione di colonna situato nella riga 'r', quindi
%occorre scambiare la riga r-esima con la riga k-esima.
%2) "STriangSup()", che risolve il sistema lineare triangolare superiore
%ignorando la triangolare inferiore.
%L'errore commesso prendendo come soluzione "x_sol" si ottiene sottraendo in
%valore assoluto "x_true" ad "x_sol", questo errore viene salvato in una
%varibile chiamata "errore".
x_sol = Sgauss_pivot(A,b);
errore=abs(x_sol-x_true);

%Punto 9.
%Calcolo la media con la formula " m=(a+b)/2 ", successivamente parte un ciclo
%for che scorre tutte le righe, quando la variabile contatore j è minore di
%3 (quindi il ciclo viene ripetuto 2 volte) con un vettore di appoggio
%salvo l'intera riga, ed effettua lo scambio tra la j-esima riga con la
%lung-j-esima riga. All'interno di questo ciclo troviamo un altro ciclo for
%che invece scorre tutte le colonne, e se trova un valore minore della
%media, allora viene sostituito con il valore della la media.
media=(L+R)/2;
lungA=length(A)+1;
for j=1:lungA-1
    if(j<3)
        Vet_App=A(lungA-j,:);
        A(lungA-j,:)= A(j,:);
        A(j,:)=Vet_App;
    end
    for l=1:lungA-1
        if(A(j,l)<media)
            A(j,l)=media;
        end
    end
end

%Punto 10.
%In questo punto definisco una vettore "x", inizialmente contenente 10 
%elementi tutti zeri, e nel ciclo for successivo assegno elemento di 
%indice "cont" di x, l'elemento della posizione "ind" del vettore X, poi 
%incremento di 2 la variabile "ind" per non far capitare elementi
%successivi.
%Dopo definisco il minimo e il massimo di x, un altro vettore
%"xx" come griglia nell'intervallo tra quel minimo e quel massimo trovati, il
%numero di punti è assegnato alla variabile "N" che in questo caso è 180.
x=zeros(1,10);
ind=1;
for cont=1:10
    x(cont)=X(ind);
    ind=ind+2;
end
mx=min(x);
Mx=max(x);
N=180;
xx = linspace(mx,Mx,180);

%Punto 11.
%Per trovare la valutazione del polinomio interpolante "pol_int" nei punti
%scelti in precedenza, utilizzo la funzione "polyval()", dova al suo 
%interno troviamo un'altra funzione, "polyfit()"; questa
%function calcola direttamente i coefficienti del polinomio, e come 
%argomenti ha bisogno: del vettore delle ascisse, del vettore delle 
%ordinate e del grado del polinomio (quindi passiamo x, y e la lunghezza di 
%x-1. Si passa la lunghezza -1 perchè il teorema dice che esiste un unico
%polinomio di grado al più n-1.
%Per trovare invece la valutazione del polinomio approssimato di grado 3, 
%basta passare come ultimo argomento 3 invece che la lunghezza-1.
%Per trovare invece la valutazione del polinomio trigonomentrico (trig) 
%ho prima di tutto creato una matrice dei coefficientri (m_c_trig) 
%composta dalla prima colonna solo di 1, dalla seconda colonna dal sin(x')
%e la terza colonna dal sin(2*x'). Ho trasposto il vettore x durante
%la crearezione della matrice in modo da trovarmi una matrice con 3 colonne
%e 10 righe. Dopo ho trovato i coefficienti della funzione trigonometrica 
%di secondo grado utilizzando la formula " matr_coef = B\y' ", così 
%facendo avremo un vettore colonna composto da i coefficienti della funzione.
%Infine ho calcolato l'approssimazione sommando il primo coefficiente, con
%il secondo coefficiente * il sin(xx) e il terzo coefficiente * il sin(2*xx).
%Purtroppo non sono riuscito a risolvere il problema di mal condizionamento 
%che da Matlab.
y = f(x);
pol_int = polyval(polyfit(x,y,length(x)-1),xx);
pol_app = polyval(polyfit(x,y,3),xx);
m_c_trig=[ones(size(x')) sin(x') sin(2*x')];
c_f_trig = m_c_trig\y';
trig = c_f_trig(1)+c_f_trig(2)*sin(xx)+c_f_trig(3)*sin(2*xx);

%Punto 12.
%In questo punto ho semplicemente mostrato graficamente con un unica
%funzione di plot, le tre interopolazioni eseguite prima, creando una nuova
%figura tramite "figure".
figure
plot(x, y,'o',xx,trig,xx,pol_app,xx,pol_int);
hold on
grid on
xlabel('Ascisse');
ylabel('Ordinate');

%Punto 13.
%Per svolgere questo punto non ho fatto altro che richiamare la function
%creata da me (media_pesata_001914(x)), in modo da definire x_mediato come
%output della function.
x_mediato= media_pesata_001914(x);

%Punto 14.
%La prima cosa è quella di generare un numero casuale reale grazie "rand()",
%successivamente assegno al vettore riga l'outuput della function creata da
%me, saxpy_001914(x,y,a).
a=rand();
sxpy = saxpy_001914(x,y,a);


%Punto 15. (prima function)
%Questa mia function prende in ingresso un solo valore, che è il vettore
%"x",e da in output un solo valore, il vettore contenente le medie pesate.
%Prima di tutto calcolo la lunghezza del vettore "x", e creo un vettore riga
%di 10 elementi contenenti tutti 0. Successivamente inizia un ciclo for che
%scorre tutta la lunghezza del vettore x. Dopo troviamo 3 if consecutivi:
%1)Controlla se si tratta del primo valore della x, se si allora calcola il
%valore con la formula (x(1)+x(2))/2.
%2)Controlla se si tratta dell'ultimo elemento della x, se si allora si
%calcola il valore con la formula (x(n-1)+x(n))/2.
%3)Controlla se il valore è compreso tra il secondo valore e il penultimo,
%se si allora calcola il valore con la formula
%(x(var_cont-1)+x(var_cont)*2+x(var_cont+1))/4;
function mx = media_pesata_001914(x)
 n=length(x);
 mx=zeros(1,10);
 for var_cont=1:n
     if(var_cont==1)
         mx(var_cont)=(x(1)+x(2))/2;
     end
     if(var_cont==n)
         mx(var_cont)=(x(n-1)+x(n))/2;
     end
     if(var_cont>1 && var_cont<n)
        mx(var_cont)=(x(var_cont-1)+x(var_cont)*2+x(var_cont+1))/4;
     end
 end
end

%Punto 16. (seconda function)
%Questa function prende in input 3 valori, il vettore "x", il vettore "y" e il
%numero casuale generato nel punto14, successivamente calcoliamo la
%lunghezza di uno dei vettori, in questo caso x. 
%Dopodichè inizia un ciclo for che scorre tutta la lunghezza del vettore e
%calcola i-esimo elemento di s con la formula: " s(i) = a*x(i) + y(i). " 
function s = saxpy_001914(x,y,a)
 lung_vett=length(x);
 s=zeros(0,10);
 for cont2=1:lung_vett
     s(cont2)=(a*x(cont2))+y(cont2);
 end
end