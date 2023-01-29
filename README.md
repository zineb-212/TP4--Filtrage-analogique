# TP4--Filtrage analogique

 ## Objectif
- Appliquer un filtre réel pour supprimer les composantes indésirables d’un signal. 

- Améliorer la qualité de filtrage en augmentant l’ordre du filtre.

### **1. Filtrage et diagramme de Bode :**
Sur le réseau électrique, un utilisateur a branché une prise CPL (Courant Porteur en Ligne), les signaux utiles sont de fréquences élevées. Le réseau électrique a cependant sa propre fréquence (50 hz). Le boiter de réception doit donc pouvoir filtrer les basses fréquences pour s'attaquer ensuite à la démodulation du signal utile.

<img width="240" alt="a"  src="https://user-images.githubusercontent.com/121026257/215289489-25b429ae-baee-4602-902f-116e58107673.PNG">

Mathématiquement, un tel filtre fournit un signal de sortie en convoluant le signal d'entrée par la réponse temporelle du filtre : *y(t) = x(t) * h(t)*
Nous souhaitons appliquer un filtre passe-haut pour supprimer la composante à 50 Hz. 

Soit notre signal d'entrée : *x(t) = sin(2.pi.f1.t) + sin(2.pi.f2.t) + sin(2.pi.f3.t)* 

Avec f1 = 500 Hz, f2 = 400 Hz et f3 = 50 Hz


- Tracage du signal
```matlab
te=1e-4;
fe=1/te;
t=0:te:5-te;
N=length(t);
f1=500;
f2=400;
f3=50;
x=sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
plot(t,x)
title('signal x(t)')
```
<img width="811" alt="1" src="https://user-images.githubusercontent.com/121026257/215290846-94e07c26-e0eb-4238-abfb-0be90216ba84.PNG">


 - Tracage du spectre.
 
 > Pour Te=0.0001
 
```matlab
 fshift = (-N/2:N/2-1)*(fe/N);
 y = fft(x);
 plot(fshift, fftshift(abs(y)/N)*2);
 title('spectre');
```
<img width="810" alt="2" src="https://user-images.githubusercontent.com/121026257/215290853-a3979305-0aa2-4387-8f2f-93d98fb133bb.PNG">

 > Pour Te=0.0005

```matlab
 %%
 te =5e-4 ;
 y = fft(x);
 plot(fshift, fftshift(abs(y)/N)*2);
 title('le spectre');
```

<img width="811" alt="3" src="https://user-images.githubusercontent.com/121026257/215290857-e1c752a2-eb2a-4bfb-af3f-843db924d8a8.PNG">

 > l'augmentation du pas d'échantillonnage rend les détails temporels du signal plus précis.

La fonction H(f) (transmittance complexe) du filtre passe haut de premier ordre est donnée par :*H(f) = (K.j.w/wc) / (1 + j. w/wc)*.

- Tracage du module de la fonction H(f).

```matlab
te = 1e-4 ;
fe = 1/te ;
t = 0:te:5 ;
N = length(t);
f = (0:N-1)*(fe/N);
K = 1 ;
wc1 = 50 ;
w = 2*pi*f ; 
H = (K*1j*w/wc1)./(1+1j*w/wc1) ;
MH = abs(H);
plot(w,MH, 'bl')
```

<img width="798" alt="4" src="https://user-images.githubusercontent.com/121026257/215294452-80d91596-5d79-4ee4-bedf-05bd24268787.PNG">

 > Le graphique obtenu représente  la façon dont le filtre affecte les différentes fréquences d'entrée. Les fréquences en dessous de la fréquence de coupure seront amplifiées, alors que les fréquences au-dessous de la fréquence de coupure seront atténuées.


- Tracage de 20.log(|H(f)|) pour différentes pulsations de coupure wc.

```matlab
te = 1e-4 ;
fe = 1/te ;
N = length(t);
f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*fe/N;
k=1;
w = 2*pi*f; 
fc1 = 30; 
fc2 =200;
fc3=300;
H1 = (k*1j*w/fc1) ./ (1 + 1j*w/fc1);
H2 = (k*1j*w/fc2) ./ (1 + 1j*w/fc2);
H3 = (k*1j*w/fc3) ./ (1 + 1j*w/fc3);
G1 = 20*log(abs(H1));
G2 = 20*log(abs(H2));
G3 = 20*log(abs(H3));
semilogx(f,G1,'g',f,G2,'r',f,G3,'b')
legend('fc = 200 rad/s', 'fc = 500 rad/s', 'fc = 1000 rad/s');
grid on 
```
<img width="801" alt="5" src="https://user-images.githubusercontent.com/121026257/215294766-ad42fe58-494f-4046-8fc4-d2c7c44dc555.PNG">


  > Nous pouvons remarquer que  plus la pulsation de coupure est petite, plus le filtre est passe haut.
- Choisix de différentes fréquences de coupure pour le domaine frequentiel 

```matlab
yt1 = y.*H1 ;
yt2 = y.*H2 ;
yt3 = y.*H3 ;
YT1 = ifft(yt1,"symmetric");
YT2 = ifft(yt2,"symmetric");
YT3 = ifft(yt3,"symmetric");
YT1_temp = fft(YT1);
YT2_temp = fft(YT2);
YT3_temp = fft(YT3);
subplot(2,2,1) 
plot(fshift,2*fftshift(abs(y))/N);
title('spectre d amplitude');
subplot(2,2,2)
plot(fshift,2*fftshift(abs(YT1_temp))/N,'r');
title('fc=200rad/s');
subplot(2,2,3)
plot(fshift,2*fftshift(abs(YT2_temp))/N,'bl');
title('fc=500rad/s');
subplot(2,2,4)
plot(fshift,2*fftshift(abs(YT3_temp))/N,'g');
title('fc=1000rad/s');
```
<img width="825" alt="6" src="https://user-images.githubusercontent.com/121026257/215295790-78b986fd-874f-4340-8f9e-e58e46c97c18.PNG">

 > A partir du code, on peut observer que ce filtre passe-haut est utilisé pour réduire les composantes de fréquences basses dans le signal en utilisant différentes fréquences de coupure . Ceci se fait par la multiplication de "y" par des filtres de fréquence de coupure differentes pour obtenir des signaux filtrés,  on applique ensuite la transformée de Fourier inverseCes sur cess ignaux  pour obtenir des signaux temporels. Les signaux temporels filtrés sont ensuite transformés à nouveau en utilisant la transformée de Fourier pour obtenir les spectres de fréquence correspondants Les résultats de ces filtrages sont ensuite tracés pour une visualisation.
 
 - Choix optimale de la pulsation de coupure wc qui vous semble optimal. 
***
```matlab
K = 1 ;
f_op=1000;
wc_op=2*pi*f_op;
w = 2*pi*f ; 
H = (K*1j*w/wc_op)./(1+1j*w/wc_op) ;
yt = y.*H ;
YT = ifft(yt,"symmetric");
YT_temp = fft(YT);
plot(fshift,2*fftshift(abs(YT_temp))/N,'r');
title('spectre filtré ');
```
<img width="809" alt="7" src="https://user-images.githubusercontent.com/121026257/215296225-a2de61e9-52bb-426f-8d99-0e3a413f9cc5.PNG">

> le filtre est bien choisi, il est observable qu'il y a une grande attenuation des freqences indesirables.


- Comparaison entre les deux signaux. 

Un filtre passe-haut ne peut jamais supprimer complètement toutes les informations indésirables car il selectionne toute une plage de fréquences . Il y aura toujours des composantes de fréquences indesirables qui passeront à travers ce filtre. La suppression de certaines fréquences peut affecter le signal.


### ** Dé-bruitage d'un signal sonore:**
 
Dans son petit studio de CROUS, un mauvais futur ingénieur a enregistré une musique en « .wav » avec un très vieux micro. Le résultat est peu concluant, un bruit strident s'est ajouté à sa musique. Heureusement son voisin, expert en traitement du signal est là pour le secourir :
*C'est un bruit très haute fréquence, il suffit de le supprimer. » dit-il sûr de lui.*
- Suppression du bruit dans ce signal

Nous devons tout dabord visualiser le signal dans le domaine frequentiel
```matlab
[y, fs] = audioread('test.wav');
%%sound(music,fs);
y = y';
N=length(y);
te = 1/fs;
t = (0:N-1)*te;
N = length(y);
f = (0:N-1)*(fs/N);
fshift = (-N/2:N/2-1)*(fs/N);
spectre_music = fft(y);
subplot(2,1,1)
plot(t, y)
title('le signal');
subplot(2,1,2)
plot(fshift,fftshift(abs(spectre_music)));
title('le spectre');
```
<img width="800" alt="8" src="https://user-images.githubusercontent.com/121026257/215296910-54d407c0-8c4d-46f7-88ae-af8e4573e9bb.PNG">

Il est evident qu'il va nous falloir un filtre passe-bas

```matlab
k = 1;
fc = 4500;
%la transmitance complexe 
h = k./(1+1j*(f/fc));
h_filter = [h(1:floor(N/2)), flip(h(1:floor(N/2)))];
y_filtr = spectre_music(1:end-1).*h_filter;
sig_filtred= ifft(y_filtr,"symmetric");
semilogx(f(1:floor(N/2)),abs( h(1:floor(N/2))),'linewidth',1.5)
```
> Nous pouvons remarquer qu'il ya eu une legere attenuation dans les haute frequences qui sont passer de plus de 20000Hz a 13000hz.

<img width="818" alt="9" src="https://user-images.githubusercontent.com/121026257/215297192-df56909f-e14e-49df-8087-e51303d6efca.PNG">

> Avant filtrage

https://user-images.githubusercontent.com/121026257/215297362-8ad710cb-90c2-414b-b015-2e5268951219.mp4


> Apres filtrage

https://user-images.githubusercontent.com/121026257/215297369-b091f402-7b65-40b2-ada2-9d64be038fa5.mp4

Nous avons créé un filtre passe-bas en utilisant la réponse impulsionnelle . La fréquence de coupure a ete fixé à 4500Hz. Ensuite nous avons utilisé une fonction de transfert complexe pour créer la réponse impulsionnelle du filtre, en un gain à 1 et un ordre de 1 pour la réponse impulsionnelle. Ensuite, nous avons multiplié le spectre du son par la réponse impulsionnelle du filtre pour obtenir le signal filtré. Nous avons ensuite appliquer la transformée de Fourier inverse pour obtenir le signal temporel filtré. 

 - L'influence du parametre K.

 Ce parametre est utilisé pour ajuster l'amplitude des fréquences conservées par le filtre. Dans notre cas nous avons k=1, ce qui signifie qu'il n'y a pas de gain en amplitude appliqué aux fréquences qui passent à travers le filtre. Alors , lorsqu'on augmentera le gain, on aura remarque une meilleur qualité du signal .

-  Amélioration de la qualité du filtrage. 

```matlab
%% Amelioration du filtre
k = 10;
fc = 4500;
%la transmitance complexe 
h = k./(1+1j*(f/fc).^100); % ordre = 100
h_filter = [h(1:floor(N/2)), flip(h(1:floor(N/2)))];
y_filtr = spectre_music(1:end-1).*h_filter;
sig_filtred= ifft(y_filtr,"symmetric");
semilogx(f(1:floor(N/2)),abs( h(1:floor(N/2))),'linewidth',1.5)
%%
plot(fshift(1:end-1),fftshift(abs(fft(sig_filtred))));
sound(sig_filtred,fs);
```

> Avant filtrage

https://user-images.githubusercontent.com/121026257/215297362-8ad710cb-90c2-414b-b015-2e5268951219.mp4


> Apres filtrage

https://user-images.githubusercontent.com/121026257/215297938-6472dd62-1649-40ca-bef9-38e01c0091d5.mp4

en augmentant l'ordre du filtre on ajoute des termes supplémentaires à l'équation qui définit le filtre. Ce qui ameliore le filtre, c'est-à-dire que les fréquences qui ne sont pas passées par le filtre seront réduites avec plus d'efficacité. Cependant, augmenter l'ordre du filtre peut également entraîner une augmentation de la distorsion du signal, car il peut y avoir une impact sur les fréquences proches de la fréquence de coupure. Dans ce cas, lorsque le paramètre K est elevé, on observe une augmentation de l'amplitude du signal final, ce qui peut rendre le son plus fort.

