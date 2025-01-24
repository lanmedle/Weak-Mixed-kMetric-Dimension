# Weak mixed k-metric dimension

## Uvod
Najina naloga je bila, da napiševa CLP program za izračun šibke mešane $k$-metrične dimenzije $wmdim_k(G)$ in določiti $\kappa''(G)$ ter $wmdim_k(G)$ za cikle, polne grafe, dvodelne polne grafe, hiperkocke in
kartezične produkte ciklov. S pomočjo tega, sva poskusila uganit bolj splošne formule za $\kappa''(G)$ in $wmdim_k(G)$.

Nato sva s pomočjo sistematičnega in stohastičnega iskanja (hill-climbing in simulated annealing) iskala grafe za katere je $wmdim_k(G)$ velik oz. majhen.

## CLP program
Kot že rečeno sva najprej morala napisati CLP program, ki bo vračal $wmdim_k(G)$. Ker mora najin program izračunati tudi razdalje med vozlišči in povezavami sva zato naprej definirala funkcijo razdalja, ki naredi prav to. Funkcija `razdalja` sprejmetri argumente: `moznost`, `vozlisce` in `G`, kjer je moznost podana kot par `(c,a)`, kjer `c` pove ali imamo povezavo ali vozlišče, `a` pa je potem ta povezava oz. vozlišče, do katerega raučnamo razdaljo. Razdalja od vozlišča do povezave `UV` je definirana kot najkrajša razdalja od vozlišča `vozlisce` do `U` oz. `V`.
```py
def razdalja(moznost, vozlisce, G):
    c, a = moznost
    if c == 'e':
        U, V = a
        return min(G.distance(U, vozlisce), G.distance(V, vozlisce))
    else:
        return G.distance(a, vozlisce)
```

CLP program za šibko mešano k-metrično dimenzijo je zelo podobne CLP programu za šibko k-metrično dimenzijo s to izjemo, da moramo pri mešani upoštevati tudi vse povezave. To sva naredila tako, da sva definirala `moznosti` kot seznam povezav in vozlišč in vsako izmed njih poimenovala - da bo funkcija razdalja vedela, ali gre za vozlišče ali povezavo (problem bi drugače nastal pri kartezičnih produktih ciklov). Pri samem programu pa sva si pomagala z vgrajeno funkcijo `MixedIntegerLinearProgram` s katero sva potem dobila `wmdim_k` od grafa `G` za nek določen `k`.

```py
def CLP_weak_mixed_k_dim(G, k):
    p = MixedIntegerLinearProgram(maximization=False)
    x = p.new_variable(binary=True)

    V = G.vertices()
    E = G.edges(labels=False)

    moznosti = [('v', v) for v in V] + [('e', e) for e in E]

    p.set_objective(sum(x[v] for v in V))
    
    for a, b in Combinations(moznosti, 2):
        p.add_constraint(
            sum(abs(razdalja(a, v, G) - razdalja(b, v, G)) * x[v] for v in V) >= k
            )

    wmdim_k = p.solve()
    mnozica_S = [v for v in V if p.get_values(x[v]) > 0.5]

    return (wmdim_k, mnozica_S)
```
Sedaj nam za prvi del preostane samo še izračun $\kappa''(G)$. Tega dobimo tako, da pogledamo do katerega `k` funkcija `CLP_weak_mixed_k_dim` ne dobi napake. Največji tak `k` je ravno $\kappa''(G)$.
```py
def kappa_2_crti(G):
    k = 1
    while True:
        try:
            CLP_weak_mixed_k_dim(G, k)
            k += 1
        except:
            return k - 1
```

## Ugotovitve
Po poganjanju kode (glej datoteko _prva_naloga.ipynb_) sva prišla do naslednjih ugotovitev.
### Cikli

$$\kappa''(G) = \Big\lfloor\frac{n}{2}\Big\rfloor$$
$$wmdim_1(G) = 3$$
$$wmdim_k(G) = \begin{cases}
    2k,& \text{če}\ n \ \text{sod}\\
    2k+1,& \text{če}\ n \ \text{lih}
\end{cases}$$

### Polni grafi
$$\kappa''(G) = 1$$
$$wmdim_1(G) = n$$

### Dvodelni polni grafi
- Za dvodelne polne grafe za katere je bil $m$ ali $n$ enak $1$, je $\kappa''(G) = 1$, $wmdim_k(G)$ je tedaj enak $m\cdot n$. Izjema je dvodelni polni graf, kjer sta $m$ in $n$ enaka $1$, tedaj je $\kappa''(G) = 1, \ wmdim_k(G) = 2$, torej $m+n$. 
- Za dvodelne polne grafe, za katere $m$ in $n$ nista enaka $1$, je $\kappa''(G) = 2$. $wmdim_k(G)$ je za $k=1$ enak $m+n-1$, če je $m$ ali $n$ enak $2$, v nasprotnem primeru pa je enak $m+n-2$. Za $k=2$ je $wmdim_k(G)$ enak $m+n$.


### Hiperkocke

$$\kappa''(G) = 2^{n-1}$$
Za $wmdim_k(G)$ žal nisva našla vzorca.

### Kartezični produkti ciklov
Pri kartezičnih produktih ciklov je za grafe s $k$-jem enakim $\kappa''(G)$ $$wmdim_k(G) = m\cdot n.$$

Če gledamo kartezični produkt ciklov $C_m$ in $C_n$, velikosti $m$ in $n$ pa je $\kappa''(G)$ definiran s pomočjo: $a = \max\{m,n\}$ in $b = \min\{m,n\}$
- če je $a$ sod: $$\kappa''(G) = \begin{cases}
\Big(\Big\lfloor \frac{b}{a} \Big\rfloor + 1\Big) \cdot a & b \ \text{lih}\\
\Big\lfloor \frac{b}{2} \Big\rfloor \cdot a& b\ \text{sod}
\end{cases}$$
- če je $a$ lih:
$$\kappa''(G) = \begin{cases}
a \cdot \Big\lfloor \frac{b}{2} \Big\rfloor - \Big\lfloor \frac{b}{2} \Big\rfloor & b \ \text{sod}\\
\Big(\Big\lfloor \frac{b}{a} \Big\rfloor + 1\Big) \cdot a& b \ \text{lih}
\end{cases}$$
- če sta $m$ in $n$ enaka: $$\kappa''(G) = \Big\lfloor \frac{m}{2} \Big\rfloor \cdot m$$




## Grafi z malim $wmdim_k(G)$
Za iskanje grafov z malim, torej $1,\, 2$ ali $3$, $wmdim_k(G)$ sva napisala kodo, ki jih bo generirala s pomočjo funkcije `nauty_geng`.
```py
def poisci_grafe_z_wmdim_k_n(od, do, n):
    for i in range(od, do + 1):
        print(f'Povezani grafi na {i} vozliscih z wmdim_k(G) = {n}:')
        for G in graphs.nauty_geng(f'{i} -c'):
            kappa_2crti = kappa_2_crti(G)
            for k in range(1, kappa_2crti + 1):
                wmdim_k, _ = CLP_weak_k_dim(G, k)
                if wmdim_k == n:
                    G.show()
```
Funkcija `poisci_grafe_z_wmdim_k_n(od, do, n)` vrača grafe (slike) od velikosti `od` do vključno velikosti `do`.  Za vsakega od njih pri vsakem `k` preveri, ali je slučajno `wmdim_k == n` in če se to zgodi, ga vrne.

### $wmdim_k(G) = 1$
Takih grafov, glede na najine rezultate ni.

### $wmdim_k(G) = 2$
Rezultat tega so poti. Spodaj so prikazani nekateri rezultati.

<table>
  <tr>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm2petvozl.png?raw=true" alt="Slika 1" width="300"/>
      <br> Slika 1: <i>Graf na petih vozliščih</i>
    </td>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm2sestvozl.png?raw=true" alt="Slika 2" width="300"/>
      <br> Slika 2: <i>Graf na šestih vozliščih</i>
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm2sedemvozl.png?raw=true" alt="Slika 3" width="300"/>
      <br> Slika 3: <i>Graf na sedmih vozliščih</i>
    </td>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm2osemvozl.png?raw=true" alt="Slika 4" width="300"/>
      <br> Slika 4: <i>Graf na osmih vozliščih</i>
    </td>
  </tr>
</table>

### $wmdim_k(G) = 3$
Grafov, ki imajo $wmdim_k(G) = 3$ je zelo veliko. Tukaj bi podala slike za grafe na petih vozliščih. Grafi na več kot petih vozliščih se nahajajo v datoteki "**druga_naloga.ipynb**"

<table>
  <tr>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm3pet1.png?raw=true" alt="Slika 5" width="300"/>
      <br> Slika 5
    </td>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm3pet2.png?raw=true" alt="Slika 6" width="300"/>
      <br> Slika 6
    </td>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm3pet3.png?raw=true" alt="Slika 7" width="300"/>
      <br> Slika 7
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm3pet4.png?raw=true" alt="Slika 8" width="300"/>
      <br> Slika 8
    </td>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm3pet5.png?raw=true" alt="Slika 9" width="300"/>
      <br> Slika 9
    </td>
    <td>
      <img src="https://github.com/lanmedle/Weak-Mixed-kMetric-Dimension/blob/main/slike/wm3pet6.png?raw=true" alt="Slika 10" width="300"/>
      <br> Slika 10
    </td>
  </tr>
</table>




## Grafi z velikim $wmdim_k(G)$
Podobno kot za manjše $wmdim_k(G)$ sva napisala tudi funkcijo za večje, torej 
```py
def poisci_grafe_z_wmdim_k_n_minus_k(od, do, k):
    for i in range(od, do + 1):
        print(f'Povezani grafi na {i} vozliscih z wmdim_k(G) = {i-k}:')
        for G in graphs.nauty_geng(f'{i} -c'):
            kappa_2crti = kappa_2_crti(G)
            for k in range(1, kappa_2crti + 1):
                wmdim_k, _ = CLP_weak_k_dim(G, k)
                if wmdim_k == i - k:
                    G.show()
```
**DODAVA SLIKE**

## Stohastično iskanje $wmdim_k(G)$ za velike grafe
```py
def generate_random_connected_graph(n, p):
    
    while True:
        G = graphs.RandomGNP(n, p)
        if G.is_connected():
            return G
```
```py
def simulated_annealing(n, k, target_wmdim, max_iterations=1000, initial_temp=10, cooling_rate=0.95):
    
    current_graph = generate_random_connected_graph(n, random.uniform(0.2, 0.8))
    current_temp = initial_temp
    
    def objective(graph):
        try:
            wmdim_k, _ = CLP_weak_mixed_k_dim(graph, k)
            return abs(wmdim_k - target_wmdim), wmdim_k
        except (ValueError, MIPSolverException):
            return float('inf'), None

    current_cost, current_wmdim = objective(current_graph)
    best_graph = current_graph
    best_cost = current_cost
    best_wmdim = current_wmdim

    for iteration in range(max_iterations):
        new_graph = current_graph.copy()
        if random.random() < 0.5:
            u, v = random.sample(range(n), 2)
            if not new_graph.has_edge(u, v):
                new_graph.add_edge(u, v)
        else:
            if new_graph.size() > n - 1:
                u, v = random.choice(new_graph.edges(labels=False))
                new_graph.delete_edge(u, v)
        
        if not new_graph.is_connected():
            continue

        new_cost, new_wmdim = objective(new_graph)
        if new_cost < current_cost or random.random() < math.exp((current_cost - new_cost) / current_temp):
            current_graph = new_graph
            current_cost = new_cost
            current_wmdim = new_wmdim
            if current_cost < best_cost:
                best_graph = current_graph
                best_cost = current_cost
                best_wmdim = current_wmdim

        current_temp *= cooling_rate
        print(f"Iteration {iteration + 1}: Current wmdim_k = {current_wmdim}, Best wmdim_k = {best_wmdim}")

        if best_cost == 0:
            break

    return best_graph, best_wmdim
```