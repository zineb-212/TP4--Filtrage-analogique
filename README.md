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


![11](https://user-images.githubusercontent.com/106840796/215270683-d5742f67-40a3-43cc-981d-9b6c80eeaf3a.PNG)
![22](https://user-images.githubusercontent.com/106840796/215271304-68528e84-9f40-43ea-b378-5fdc63b23a12.PNG)
![33](https://user-images.githubusercontent.com/106840796/215271302-49052c7f-70fa-44b8-9bc9-e8a6ac9b14d8.PNG)
***
 ### **Explication :**
 ###### on a dabord créé un filtre passe-bas en utilisant la réponse impulsionnelle . La fréquence de coupure (fc) est définie comme 4500Hz. et puis , on a utilisé une fonction de transfert complexe pour créer la réponse impulsionnelle du filtre, en utilisant un paramètre k=1 pour le gain et un ordre de 1 pour la réponse impulsionnelle. Ensuite, on a multiplié le spectre de la musique par la réponse impulsionnelle du filtre pour obtenir le signal filtré. et on a utilisé ensuite la transformée de Fourier inverse pour obtenir le signal temporel filtré. On observe qu on n a pas encore pu supprimer toutes les frequences indesirables .

***
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ [ (Revenir au sommaire) ](#retour)
***
#### $~~~~~~$ **2. Mettez-la en oeuvre. Quelle influence à le paramètre K du filtre que vous avez utilisé ?** 
***
```matlab
% qst 2
%on a augmenté la valeur de K :
k = 10;
fc = 4500;
%la transmitance complexe 
h = k./(1+1j*(f/fc));

h_filter = [h(1:floor(N/2)), flip(h(1:floor(N/2)))];
y_filtr = spectre_music(1:end-1).*h_filter;
sig_filtred= ifft(y_filtr,"symmetric");
semilogx(f(1:floor(N/2)),abs( h(1:floor(N/2))),'linewidth',1.5)
title('le gain =10 :');
%%
plot(fshift(1:end-1),fftshift(abs(fft(sig_filtred))));
%sound(sig_filtred,fs);
title('le spectre du signal filtré:');
```
***
![gain10](https://user-images.githubusercontent.com/106840796/215271677-76bc84fe-27ca-4707-96ea-1cd8a55739ea.PNG)
![10](https://user-images.githubusercontent.com/106840796/215271663-3ed26c84-25b0-4e7b-b6b3-faba14dfd97c.PNG)

 ### **Explication :**
 ###### Le paramètre K du filtre est le gain en amplitude à la fréquence de coupure (fc) du filtre. Il est utilisé pour augmenter ou diminuer l'amplitude des fréquences qui passent à travers le filtre. Plus le paramètre K est élevé, plus l'amplitude des fréquences autour de la fréquence de coupure sera élevée, et inversement. Il est utilisé pour ajuster l'amplitude des fréquences conservées par le filtre. Dans le code spécifié, le paramètre K est égal à 1, ce qui signifie qu'il n'y a pas de gain en amplitude appliqué aux fréquences qui passent à travers le filtre. Alors , lorsqu on augmenté le gain a 10 , cad amplification , on remarque un progres du qualité du signal .
***
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ [ (Revenir au sommaire) ](#retour)
***
#### $~~~~~~$ **3. Quelles remarques pouvez-vous faire notamment sur la sonorité du signal final.** 
***
 ###### le paramètre K du filtre correspond à l'amplitude de sortie du filtre.En effet , lorsqu on a augmenter  K a 10, l'amplitude des fréquences qui passent à travers le filtre est élevée, c est l amplification .

***
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$ [ (Revenir au sommaire) ](#retour)
***
#### $~~~~~~$ **4. Améliorer la qualité de filtrage en augmentant l’ordre du filtre.** 
***
```matlab
%   qst 4
k = 10;
fc = 4500;
%la transmitance complexe 
h = k./(1+1j*(f/fc).^100);

h_filter = [h(1:floor(N/2)), flip(h(1:floor(N/2)))];
y_filtr = spectre_music(1:end-1).*h_filter;
sig_filtred= ifft(y_filtr,"symmetric");
semilogx(f(1:floor(N/2)),abs( h(1:floor(N/2))),'linewidth',1.5)
title('le gain =10  et ordre=100 :');
%%
plot(fshift(1:end-1),fftshift(abs(fft(sig_filtred))));
%sound(sig_filtred,fs);
title('le spectre du signal filtré pour ordre=100:');
```
***
![ok](https://user-images.githubusercontent.com/106840796/215272155-3582dc8d-117e-4e07-aab8-2a949748da72.PNG)
![yes](https://user-images.githubusercontent.com/106840796/215272157-96027833-0d68-4fc2-bf95-285007e32016.PNG)

 ### **Explication :**
 ###### Lorsqu'on augmente l'ordre du filtre, cela signifie qu'on ajoute des termes supplémentaires à l'équation qui définit le filtre. Cela peut entraîner une amélioration de la réjection des fréquences indésirables, c'est-à-dire que les fréquences qui ne sont pas passées par le filtre seront réduites avec plus d'efficacité. Cependant, augmenter l'ordre du filtre peut également entraîner une augmentation de la distorsion du signal, car il peut y avoir une réaction plus importante aux fréquences proches de la fréquence de coupure. Dans ce cas particulier, lorsque le paramètre K est augmenté, on observe une augmentation de l'amplitude du signal final, ce qui peut rendre le son plus fort mais avec un risque d'écrêtage (saturation) et de distorsion.Ainsi , on est plus proche du signal ideal , cad il n y a pas de perte d infos , et le bruit sera supprime ,on n aura pas la bande de transition  .
***

