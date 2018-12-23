# Per implementare viscosità
## Perchè implementare viscosità
L'idea di implementare la viscosità nasce dal fatto che quando ho poca dens. di massa appena fuori dall'uscita del capillare, la simulazione si ferma perchè il calcolo della conduzione termica diventa instabile (secondo me a causa delle perturbazioni dovuti all'avvezione, combinato col fatto che k/cv sia molto elevato). Con un po' di viscosità forse si smorza un po' il gradiente di densità, perchè la viscosità dovrebbe smussare i gradienti e le oscillazioni di velocità e far avere un comportamento migliore alla fisica attorno agli shock. Oltretutto succede anche spesso che l'energia interna diventa negativa perchè l'energia cinetica è molto più grande della en.interna. Mentre con un po' di viscosità si potrebbe convertire un po' di en.cinetica in interna.
## Come
Si può sfruttare l'implementazione STS già presente in PLUTO, basta assicurarsi che sia compatibile con l'ADI.
Penso sia compatibile, solo che vanno aggiustati alcuni dettagli.
+ Se attivo la viscosità, dato che uso le coordinate cilindriche, devo stare attento a modificare correttamente il calcolo del rhs per togliere i termini sorgenti dovuti alla viscosità in coordinate curvilinee nel caso in cui scelgo di tenere bloccato il fluido; infatti sulla guida a pag 72 c'è scritto _"In curvilinear geometries, additional geometrical source terms coming from the tensor’s divergence are added to the right hand side of the equations."_.

Come implementare all'atto pratico:
+ Devo capire bene quali sono le equazioni che vengono risolte e come si normalizzano, perchè poi devo normalizzare la viscisità correnttamente. Comunque la viscisità deve avere dimensioni $\rho \times \mathrm{length}^2/\mathrm{time}$, vedi pag 73 della guida. Ricorda che per gas monoatomici $\nu_2=0$ (vedi la def. di $\nu_2$ sulla guida).
+ Trovare e implementare un caso test per verificare che la viscosità funzioni bene.
+ Prima provare con una viscosità che sia realistica per l'idrogeno ma non accurata, se no rischio di perdere tempo con il coefficiente di visc. quando magari il codice è instabile o il limite sul dt per avanzamento STS è proibitivo. 
