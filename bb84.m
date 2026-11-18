% --- BB84 szimuláció: Alice -> zajos, csillapító csatorna -> Bob (detektorokkal, afterpulse) ---
% Sűrűségmátrixos leírás, ideális egyfoton-forrás, a csatornában bit-flip és
% depolarizáló zaj + csillapítás, Bobnál nem ideális, afterpulse-os detektorok.
%
%   Alice:    véletlen bit + bázis -> rho = |psi><psi| (a 4 állapot egyike)
%   Csatorna: Kraus-operátorok, rho -> sum_k K_k rho K_k'
%             depolarizáció: rho -> (1-p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
%             bit-flip:      rho -> (1-p) rho + p X rho X
%             csillapítás:   rho -> T rho
%   Bob mérés: POVM, E_k = eta_k * Pi_k  (detektor-hatásfok)
%             P(D0) = Tr(E_0 rho_B), P(D1) = Tr(E_1 rho_B), a maradék: nincs jelklikk
%   Detektorok (a qkdsim_v1 alapján): jitter, sötétzaj, kapuzás, holtidő,
%             dupla klikk -> eldobva
%   Afterpulse: egy lavina p_ap eséllyel csapdába ejt egy töltéshordozót, amely
%             Exp(tau_ap) idő múlva szabadul fel. Ha ez a holtidő után, egy kapun
%             belül történik, ugyanaz a detektor hamis klikket ad, ami újabb
%             afterpulse-t indíthat.
%   Sifting:  csak a detektált ÉS egyező bázisú pulzusok maradnak
%   Utófeldolgozás

clear; clc; rng('shuffle');

%% Paraméterek
f_rep        = 1e6;     % órajel [Hz]
T_sim        = 10;       % szimulált idő [s]
chunk_pulses = 1e6;     % egy blokk mérete (memória/sebesség tradeoff)

dt       = 1 / f_rep;   % pulzusperiódus [s]
N_total  = round(T_sim * f_rep);
N_chunks = ceil(N_total / chunk_pulses);

%% Csatorna — paraméter
L             = 20;     % [km]
alpha         = 0.2;    % [dB/km]
transmittance = 10^(-(alpha*L)/10);

p_bf  = 0.02;           % bit-flip valószínűsége
p_dep = 0.03;           % depolarizáció valószínűsége

%% Bob — detektorok — paraméter (D0 = bit 0, D1 = bit 1)
eta_D       = [0.10, 0.10];     % detektor-hatásfok, beérkező fotonok hányad része kattan
lambda_dark = [200, 200];       % sötétszámlálási ráta [Hz], hamis klikkek, ha nincs foton
tau_dead    = 50e-9;            % holtidő [s] (nem bénuló), a detektor egy klikk után tau_dead ideig nem képes új klikket adni, holtidő alatt érkező klikk elveszik
T_gate      = 200e-12;          % detekciós kapu szélessége [s], a pulzus körüli időablak, amikor a detektor képes klikket adni
sigma_jit   = 30e-12;           % időbeli jitter (1 sigma) [s], a klikk időpontja a pulzus közepétől eltérhet
p_ap        = [0.05, 0.05];     % P(egy lavina csapdába ejt egy töltéshordozót)
tau_ap      = 100e-9;           % a csapdázott töltéshordozó élettartama [s]

assert(T_gate <= dt, 'A kapu (T_gate) nem lehet hosszabb a pulzusperiódusnál.');

%% Utófeldolgozás — paraméter
v_pe    = 0.5;          % a paraméterbecslésre feláldozott hányad (előadás: BB84-nél a fele)
Q_tol   = 0.11;         % QBER-küszöb: e fölött a protokoll leáll
eps_pe  = 1e-10;        % P(a maradék bitek QBER-je meghaladja a becsült felső korlátot)
n_frame = 1e4;          % a Cascade egy frame-jének hossza [bit]
n_pass  = 4;            % Cascade körök száma
t_conf  = 64;           % a megerősítő hash hossza [bit]

%% Állapotok (Alice)
ket0 = [1; 0];
ket1 = [0; 1];
ketP = (ket0 + ket1) / sqrt(2);
ketM = (ket0 - ket1) / sqrt(2);
kets  = [ket0, ket1, ketP, ketM];   % állapotindex: s = 2*bázis + bit + 1, matlab indexelés miatt
nevek = {'|0>', '|1>', '|+>', '|->'};


% Sűrűségmátrixok: rho(:,:,s) = |psi_j><psi_j|, j = 1..4
rho = zeros(2, 2, 4);
for j = 1:4
    rho(:,:,j) = kets(:,j) * kets(:,j)';
end

%% Mérés (Bob)
% Projektorok: Pi(:,:,kimenet,bázis); kimenet 1 = D0 (bit 0), 2 = D1 (bit 1)
%                                     bázis   1 = Z,          2 = X
%                                     Z: |0><0|, |1><1|; X: |+><+|, |-><-|
Pi = zeros(2, 2, 2, 2);
Pi(:,:,1,1) = rho(:,:,1);  Pi(:,:,2,1) = rho(:,:,2);
Pi(:,:,1,2) = rho(:,:,3);  Pi(:,:,2,2) = rho(:,:,4);

% Nem ideális detektorok: POVM-elemek E(:,:,kimenet,bázis)
%   E_k = eta_k * Pi_k                    detektor-hatásfok
% Pi_0 + Pi_1 = I, de E_0 + E_1 < I: a maradék I - E_0 - E_1 a "nincs jelklikk"
% A csillapítás nem itt van, hanem a csatornában (K = sqrt(T) * I), így
%   P(D_k) = Tr(E_k * rho_B) = T * eta_k * Tr(Pi_k * rho)
E = zeros(2, 2, 2, 2);
for b = 1:2
    E(:,:,1,b) = eta_D(1) * Pi(:,:,1,b);
    E(:,:,2,b) = eta_D(2) * Pi(:,:,2,b);
end

%% Csatorna (Kraus-operátorok)
% Csillapítás: a foton T valószínűséggel jut át, a polarizációja nem változik
% sum_k K_k' K_k = T * I < I
% Tr(rho_B) = T a megérkezés valószínűsége, a hiányzó 1 - T az elveszett foton
% (A detektor-hatásfok a POVM-ben van: E_k = eta_k * Pi_k.)
%             depolarizáció: rho -> (1-p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
%             bit-flip:      rho -> (1-p) rho + p X rho X
%             csillapítás:   rho -> T rho
I2 = eye(2);
X  = [0 1; 1 0];                    % Pauli-mátrixok
Y  = [0 -1i; 1i 0];
Z  = [1 0; 0 -1];

% Depolarizáció:  rho -> (1-p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
K_dep  = {sqrt(1 - p_dep) * I2, sqrt(p_dep/3) * X, sqrt(p_dep/3) * Y, sqrt(p_dep/3) * Z};
% Bit-flip:       rho -> (1-p) rho + p X rho X
K_bf   = {sqrt(1 - p_bf) * I2, sqrt(p_bf) * X};
% Csillapítás:    rho -> T rho   (nyomcsökkentő: Tr(rho_B) = T a megérkezés valószínűsége)
K_loss = {sqrt(transmittance) * I2};

% Egymás utáni csatornák: a Kraus-operátorok minden párosításának szorzata,
% itt 4 x 2 x 1 = 8 operátor. A jobb szélső hat először (depolarizáció ->
% bit-flip -> csillapítás); ezek a csatornák felcserélhetők
K = {};
for a1 = 1:numel(K_loss)
    for a2 = 1:numel(K_bf)
        for a3 = 1:numel(K_dep)
            K{end+1} = K_loss{a1} * K_bf{a2} * K_dep{a3}; %#ok<SAGROW>
        end
    end
end

% a komment a SAGROW-val azt csinálja, hogy a MATLAB nem figyelmeztet, hogy a K tömb mérete változik a ciklusban, és emiatt lassabb lehet a kód

% Ellenőrzés: sum_k K_k' K_k = T * I  (a zaj nyomtartó, csak a csillapítás csökkenti)
% ha sajátérték > 1, akkor a csatorna nem nyomtartó, a QBER > 50% lesz
KK = zeros(2);
for k = 1:numel(K)
    KK = KK + K{k}' * K{k};
end

rho_B = zeros(2, 2, 4);
for j = 1:4
    for k = 1:numel(K)
        rho_B(:,:,j) = rho_B(:,:,j) + K{k} * rho(:,:,j) * K{k}';
    end
end

%% Táblázat: P0(s, b) = Tr(E_0^b * rho_B^s),  P1(s, b) = Tr(E_1^b * rho_B^s)
% A csatorna és Bob mérése minden pulzusra ugyanaz, ezért a 4 állapot x 2 bázis
% valószínűségét elég egyszer kiszámolni, a ciklusban csak kiolvassuk
% Ez a jelklikk valószínűsége; a sötétzaj és a holtidő a detektor-részben jön
% P0 + P1 = T < 1, ezért a D1 valószínűségét is külön kell tárolni
P0 = zeros(4, 2);
P1 = zeros(4, 2);
for j = 1:4
    for b = 1:2
        P0(j,b) = real(trace(E(:,:,1,b) * rho_B(:,:,j))); % számitásnál kerekítési hiba miatt kis képzetes rész lehet, ezért real() függvény
        P1(j,b) = real(trace(E(:,:,2,b) * rho_B(:,:,j)));
    end
end

fprintf('Csatorna: L = %g km, alpha = %g dB/km -> transmittance = %.3g (%.1f dB)\n', ...
    L, alpha, transmittance, alpha*L);
fprintf('Zaj: bit-flip p = %g, depolarizáció p = %g | %d Kraus-operátor, sum K''K = %.4g * I\n', ...
    p_bf, p_dep, numel(K), real(KK(1,1)));
fprintf('Detektorok: eta = [%g %g], sötét = [%g %g] Hz, holtidő = %g ns, kapu = %g ps, jitter = %g ps\n', ...
    eta_D, lambda_dark, tau_dead*1e9, T_gate*1e12, sigma_jit*1e12);
fprintf('Afterpulse: p_ap = [%g %g], tau_ap = %g ns\n', p_ap, tau_ap*1e9);
fprintf('\nJelklikk-valószínűségek pulzusonként:\n');
fprintf('Állapot | Tr(rho_B)  | Z: P(D0)   Z: P(D1)   X: P(D0)   X: P(D1)\n');
for j = 1:4
    fprintf('  %s   | %.3g  | %.3g  %.3g  %.3g  %.3g\n', nevek{j}, ...
        real(trace(rho_B(:,:,j))), P0(j,1), P1(j,1), P0(j,2), P1(j,2));
end

%% Akkumulátorok
% A sifting és QBER számításához
% tot_detected: összes detektált bit
% tot_sifted: összes siftelt bit
% tot_errors: összes siftelt bitből hány hibás
% tot_n0: összes siftelt bitből hány 0
tot_detected = 0; tot_sifted = 0; tot_errors = 0; tot_n0 = 0;
tot_sifted_b = [0 0];               % siftelt bitek bázisonként (Z, X)
tot_errors_b = [0 0];               % hibák bázisonként (Z, X)
tot_double   = 0;                   % dupla klikkek (eldobva)
tot_src      = zeros(2, 3);         % detektor x forrás (1 = jel, 2 = sötét, 3 = afterpulse)
last_click   = [-inf, -inf];        % holtidő-állapot: utolsó klikk abszolút ideje
ap_pending   = {zeros(1,0), zeros(1,0)};   % függő felszabadulások abszolút ideje detektoronként
key_a        = cell(1, N_chunks);   % Alice és Bob siftelt kulcsa chunkonként
key_b        = cell(1, N_chunks);

fprintf('\nIndulás: %d pulzus (%g s @ %g GHz), %d chunk\n', ...
    N_total, T_sim, f_rep/1e9, N_chunks);
tic;
%% Fő ciklus
for c = 1:N_chunks
    % ---- Chunk mérete -----
    % Az utolsó chunk lehet rövidebb, ha N_total nem osztható chunk_pulses-szal
    % Egy blokkopn belül minen vektorosan fut
    N       = min(chunk_pulses, N_total - (c-1)*chunk_pulses);
    t0      = (c-1) * chunk_pulses * dt;    % a chunk abszolút kezdőideje
    T_chunk = N * dt;

    % ----- Alice -----
    % Alice véletlen bitjei és bázisai
    alice_bits  = randi([0 1], 1, N);
    alice_bases = randi([0 1], 1, N);           % 0 = Z, 1 = X
    j_idx       = 2*alice_bases + alice_bits + 1;

    % ----- Bob: bázis + jelklikk -----
    bob_bases = randi([0 1], 1, N);
    index = j_idx + 4*bob_bases;                  % lineáris index: = (j_idx, bob_bases+1)
    p0  = P0(index);                              % p0 = P0(j_idx + 4*bob_bases)
    p1  = P1(index);                                  % p1 = P1(j_idx + 4*bob_bases)

    % Egy véletlen szám, három kimenet:
    %   [0, p0)        -> jelklikk D0-n
    %   [p0, p0+p1)    -> jelklikk D1-en
    %   [p0+p1, 1)     -> nincs jelklikk (elveszett a csatornában / a detektorban)
    u = rand(1, N);
    sig_idx = {find(u < p0), find(u >= p0 & u < p0 + p1)};   
    % pulzusindexek D0-ra, D1-re. A detektor-hatások miatt ezek egy része elveszhet, és sötétklikkek is hozzáadódhatnak

    % ----- Detektorok: jitter, sötétzaj, kapuzás, holtidő, afterpulse -----
    pulse_clicks = false(2, N);                 % kattant-e D0/D1 a k-adik pulzus kapujában
    for d = 1:2
        % jel: a k-adik pulzus a periódus közepén, (k - 1/2)*dt-kor érkezik, + jitter
        t_sig  = (sig_idx{d} - 0.5) * dt + sigma_jit * randn(1, numel(sig_idx{d}));
        % sötétzaj: Poisson-folyamat a teljes chunkra (a kapun kívül is):
        % a darabszám Poisson-eloszlású, az időpontok egyenletesen szórnak
        % itt kell a Statistics and Machine Learning Toolbox, a poissrnd() függvény miatt
        N_dark = poissrnd(lambda_dark(d) * T_chunk);
        t_dark = T_chunk * rand(1, N_dark);

        % a jel és sötétzaj időpontjait egyesítjük, és sorba rendezzük
        % a forrás típusát is tároljuk (1 = jel, 2 = sötét)
        t   = [t_sig, t_dark];
        src = [ones(1, numel(t_sig)), 2*ones(1, numel(t_dark))];

        % kapuzás: csak a pulzusidő körüli T_gate ablakba eső esemény számít
        k = round(t/dt + 0.5);                  % a legközelebbi pulzus indexe
        in_gate = k >= 1 & k <= N & abs(t - (k - 0.5)*dt) <= T_gate/2; % az esemény a k-adik pulzus ablakában van-e
        t = t(in_gate);  src = src(in_gate);  k = k(in_gate); % csak a kapun belüli események maradnak

        % sorba rendezzük az eseményeket idő szerint, hogy a holtidőt alkalmazhassuk
        [t, ord] = sort(t);
        src = src(ord);  k = k(ord);

        % holtidő + afterpulse: egy klikk után tau_dead ideig a detektor vak.
        % előző klikk ideje: last_click(d), a chunkon belüli idő: last_click(d) - t0
        % Időrendben kell haladni, mert minden elfogadott klikk új afterpulse-t 
        % indíthat, és az afterpulse is klikk: holtidőt és újabb afterpulse-t okoz
        % A holtidő-állapot és a függő felszabadulások a chunkok közt is megmaradnak
        n      = numel(t);                      % összes jel/sötét esemény a chunkban
        keep   = false(1, n);                   % melyik jel/sötét esemény marad meg (nem esik holtidőbe)
        t_ap   = zeros(1, 0);                   % elfogadott afterpulse-ok ideje
        last_t = last_click(d) - t0;            % abszolút -> chunkon belüli idő
        pend   = ap_pending{d} - t0;            % függő felszabadulások, növekvő sorrendben
        i = 1;

        % addig fut, amíg van jel/sötét esemény, vagy van függő afterpulse, amely a chunkban esik
        while i <= n || (~isempty(pend) && pend(1) <= T_chunk)
            % a következő esemény: a jel/sötét lista vagy a függő afterpulse, amelyik korábbi
            is_ap = ~isempty(pend) && (i > n || pend(1) < t(i));

            % ha van függő afterpulse, és vagy elfogyott a t lista (i > n), vagy a függő afterpulse korábbi, mint a következő jel/sötét esemény.
            % akkor a következő esemény az afterpulse, különben a következő jel/sötét esemény
            if is_ap
                t_e = pend(1);
                pend(1) = [];
            else
                t_e = t(i);
            end
            % t_e a következő esemény ideje (chunkon belül)



            if t_e - last_t >= tau_dead % tau_dead ideig nem észlelhet klikket
                last_t = t_e;
                if is_ap
                    t_ap(end+1) = t_e; %#ok<SAGROW>
                else
                    keep(i) = true; % a klikket megtartjuk, ha a holtidőn kívül van
                end
                % A lavina p_ap eséllyel csapdába ejt egy töltéshordozót. Ez csak
                % akkor okoz klikket, ha a holtidő után, egy ablakon belül szabadul fel.
                if rand < p_ap(d)
                    delay = -tau_ap * log(rand);        % Exp(tau_ap) eloszlású
                    t_rel = t_e + delay;
                    if delay >= tau_dead && abs(mod(t_rel, dt) - dt/2) <= T_gate/2
                        pend = sort([pend, t_rel]); % a függő afterpulse-ok listájába kerül, növekvő sorrendben
                    end
                end
            end

            % A holtidő miatt a következő jel/sötét esemény indexe csak akkor nő, ha a következő esemény nem afterpulse volt.
            % Ha az afterpulse volt a következő esemény, akkor i nem változik, és a while ciklusban újra ellenőrizzük a következő jel/sötét eseményt
            if ~is_ap
                i = i + 1;
            end
        end
        % A chunk végén a holtidő-állapot és a függő afterpulse-ok megmaradnak a következő chunkra
        last_click(d) = t0 + last_t;
        ap_pending{d} = pend + t0;

        % pulzusok kattantak-e a detektoron (jel vagy sötét), és a forrás típusát is számoljuk
        % tot_src(d,1) = jelkikkek száma, tot_src(d,2) = sötétklikkek száma, tot_src(d,3) = afterpulse-ok száma
        pulse_clicks(d, [k(keep), round(t_ap/dt + 0.5)]) = true;
        tot_src(d,:) = tot_src(d,:) + [nnz(src(keep) == 1), nnz(src(keep) == 2), numel(t_ap)];
    end

    % ----- Bob bitje -----
    %   csak D0 -> 0 | csak D1 -> 1 | egyik sem -> -1
    %   mindkettő (dupla klikk) -> a bitet eldobjuk (-1 marad), de számon tartjuk
    c0 = pulse_clicks(1,:);
    c1 = pulse_clicks(2,:);
    bob_bits = -ones(1, N);
    bob_bits( c0 & ~c1) = 0;
    bob_bits(~c0 &  c1) = 1;
    both = c0 & c1;
    detected = (bob_bits >= 0);                 % csak az egy klikkes pulzusok

    % ----- Sifting -----
    % Bob nyilvánosan bejelenti, hol volt klikk és milyen bázisban mért;
    % csak a detektált ÉS egyező bázisú pulzusok maradnak.
    sift = detected & (alice_bases == bob_bases);
    tot_detected    = tot_detected + nnz(detected);
    tot_double      = tot_double + nnz(both);
    tot_sifted   = tot_sifted + nnz(sift);
    tot_errors   = tot_errors + nnz(alice_bits(sift) ~= bob_bits(sift));
    tot_n0       = tot_n0 + nnz(bob_bits(sift) == 0);
    key_a{c}     = alice_bits(sift) == 1;
    key_b{c}     = bob_bits(sift) == 1;

    % bázisonként is: a bit-flip csak a Z bázisban okoz hibát
    for b = 1:2
        sb = sift & (bob_bases == b-1);
        tot_sifted_b(b) = tot_sifted_b(b) + nnz(sb);
        tot_errors_b(b) = tot_errors_b(b) + nnz(alice_bits(sb) ~= bob_bits(sb));
    end

    if mod(c, max(1, round(N_chunks/10))) == 0 || c == N_chunks
        fprintf('  chunk %d/%d (%.0f%%) | eltelt %.1f s\n', c, N_chunks, 100*c/N_chunks, toc);
    end
end

%% Eredmények
QBER   = tot_errors / max(tot_sifted, 1);
QBER_b = tot_errors_b ./ max(tot_sifted_b, 1);

% Várt értékek holtidő nélkül (közelítés: átlagos hatásfok)
p_sig  = transmittance * mean(eta_D);            % jelklikk / pulzus
p_dark = lambda_dark * T_gate;                   % sötétklikk / pulzus, detektoronként
p_prim = p_sig + 2*mean(p_dark);                 % elsődleges (jel + sötét) klikk / pulzus
R_tau  = f_rep * (p_prim/2) * tau_dead;          % klikkráta * holtidő, detektoronként

% Afterpulse / klikk: a felszabadulás a holtidő után, egy kapuba esik
% (a klikket a kapu közepére véve; telítésnél a detektor gyakran vak, ezért a mért
% érték kisebb). A láncokkal együtt egy elsődleges klikkre p/(1-p) afterpulse jut.
n_g   = floor(tau_dead/dt) : ceil((tau_dead + 30*tau_ap)/dt);    % a klikk utáni kapuk
lo    = max(n_g*dt - T_gate/2, tau_dead);
hi    = max(n_g*dt + T_gate/2, tau_dead);
p_eff = mean(p_ap) * sum(exp(-lo/tau_ap) - exp(-hi/tau_ap));
a_v   = p_eff / (1 - p_eff);                                     % várt
a_m   = sum(tot_src(:,3)) / max(sum(tot_src(:,1:2), 'all'), 1);  % mért

% a csatorna-zaj hibája egyező bázisnál:
q   = 2*p_dep/3;                                 % depolarizáció: X és Y fordít Z-ben, Y és Z fordít X-ben
e_b = [p_bf*(1-q) + (1-p_bf)*q, q];              % Z: a bit-flip VAGY a depolarizáció fordít (mindkettő = nincs hiba)
                                                 % X: a bit-flip nem hat (|+>, |-> az X sajátállapotai)
% Az afterpulse egy másik pulzusnál kattan, ezért a bitje véletlen: a fele hibás.
QBER_v_b = (e_b*p_sig + mean(p_dark) + a_m/2*p_prim) / ((1 + a_m)*p_prim);
QBER_v   = mean(QBER_v_b);                       % a két bázis egyforma valószínű

fprintf('\n--- BB84 (zajos csatorna + detektorok + afterpulse) eredmények ---\n');
fprintf('Szimulált idő: %g s | futásidő: %.1f s\n', N_total/f_rep, toc);
fprintf('Elküldött pulzusok: %d\n', N_total);
for d = 1:2
    fprintf('  D%d klikkek: %d  [jel = %d, sötét = %d, afterpulse = %d]\n', d-1, sum(tot_src(d,:)), ...
        tot_src(d,1), tot_src(d,2), tot_src(d,3));
end
fprintf('Afterpulse / elsődleges klikk: %.4f%%  (várt telítés nélkül: %.4f%%)\n', 100*a_m, 100*a_v);
fprintf('Detektált (egy klikkes) pulzusok: %d  (várt holtidő nélkül: %.0f)\n', ...
    tot_detected, N_total*p_prim*(1 + a_v));
fprintf('Dupla klikkek (eldobva): %d\n', tot_double);
fprintf('  R * tau_dead = %.3g detektoronként (ha nem << 1, a holtidő sok klikket elvesz)\n', R_tau);
fprintf('Siftelt kulcs: %d bit | siftelt kulcsráta: %.4g bit/s\n', ...
    tot_sifted, tot_sifted/(N_total/f_rep));
fprintf('QBER:   %.4f%%  (várt a mért afterpulse-aránnyal: %.4f%%)\n', 100*QBER, 100*QBER_v);
fprintf('  Z-ben: %.4f%%  (várt: %.4f%%)\n', 100*QBER_b(1), 100*QBER_v_b(1));
fprintf('  X-ben: %.4f%%  (várt: %.4f%%)\n', 100*QBER_b(2), 100*QBER_v_b(2));
fprintf('Bias: P(0) = %.6f\n', tot_n0 / max(tot_sifted, 1));

%% Utófeldolgozás
% A klasszikus csatornát hitelesítettnek tekintjük
fprintf('\n--- Utófeldolgozás ---\n');
tic;
key_a  = [key_a{:}];
key_b  = [key_b{:}];
n_sift = numel(key_a);

% ----- 1. Paraméterbecslés -----
% A siftelt bitek véletlen részét nyilvánosan összevetik, majd eldobják.
k_pe  = round(v_pe * n_sift); % a mintába kerülő bitek száma
is_pe = false(1, n_sift); % a mintába kerülő bitek indexe
is_pe(randperm(n_sift, k_pe)) = true; % a mintába kerülő bitek indexének a véletlen kiválasztása
Q_k   = nnz(key_a(is_pe) ~= key_b(is_pe)) / max(k_pe, 1); % a mintában mért QBER
ka    = key_a(~is_pe);  % a maradék siftelt bitek (Alice)
kb    = key_b(~is_pe);  % a maradék siftelt bitek (Bob)
n_pp  = numel(ka);      % utófeldolgozásra maradó bitek száma
if k_pe == 0 || n_pp == 0 % túl kevés a siftelt bit, nincs kulcs
    fprintf('Túl kevés a siftelt bit (%d), nincs kulcs.\n', n_sift);
    return;
end

% Serfling-egyenlőtlenség: a maradék n bit a teljes N = n + k bit véletlen részhalmaza,
% ezért P(Q_n >= Q_k + beta) <= exp(-2 n k^2 beta^2 / (N (k+1))). Ebből beta:
beta = sqrt(n_sift * (k_pe + 1) * log(1/eps_pe) / (2 * n_pp * k_pe^2)); % a maradék QBER-jének felső korlátja a mintából becsülve
Q_up  = min(Q_k + beta, 1); % kulcsbitek QBER-jének felső korlátja a mintából becsülve

fprintf('1. Paraméterbecslés: %d siftelt bitből %d mintabit, %d marad\n', n_sift, k_pe, n_pp);
fprintf('   Minta QBER: %.4f%% | felső korlát a maradékra (Serfling): %.4f%%\n', 100*Q_k, 100*Q_up);
if Q_k > Q_tol
    fprintf('   A minta QBER a %.2f%%-os küszöb fölött van: a protokoll leáll, nincs kulcs.\n', 100*Q_tol);
    return;
end

% ----- 2-3. Hibajavítás (Cascade) és megerősítés, frame-enként -----
n_frames   = ceil(n_pp / n_frame);  % a frame-ek száma a hibajavításnál
kept       = false(1, n_pp);        % a megerősítésen átment bitek
leak_kept  = 0;                     % a megtartott frame-ek nyilvános paritásbitjei
n_kept     = 0;                     % a megtartott frame-ek száma
n_err_pre  = nnz(ka ~= kb);         % a hibák száma a Cascade előtt
for f = 1:n_frames
    idx = (f-1)*n_frame + 1 : min(f*n_frame, n_pp); % az f-edik frame indexei
    fa  = ka(idx);                  % a frame Alice-nál
    fb  = kb(idx);                  % a frame Bobnál, ezt javítjuk
    n_f = numel(idx);               % a frame hossza
    if n_f == 0
        % TODO: a frame lehet 0 hosszú?
        continue;
    end

    % --- Cascade: Bob kijavítja fb-t Alice fa-jára ---
    % Egy kör: véletlen keverés, blokkokra osztás, a páratlan blokkokban bináris
    % kereséssel a hiba megkeresése. Cascade-hatás: a kijavított bit a korábbi körök
    % blokkjainak paritását is megfordítja, így ott újabb hibák kereshetők meg.
    % leak = a nyilvánosan elküldött paritásbitek száma. Az fa-t a kód csak a
    % paritásokhoz olvassa, ezek mindegyike beleszámít a leak-be.
    perm  = cell(1, n_pass);        % perm{p}: a p-edik kör keverése
    pos   = cell(1, n_pass);        % pos{p}(i): az i-edik bit helye a keverés után
    odd   = cell(1, n_pass);        % odd{p}(j): a j-edik blokk paritása eltér-e
    k_blk = zeros(1, n_pass);       % blokkméret körönként: 0.73/Q_k, utána duplázódik
    leak  = 0;

    % egy kör: a blokkméret a mintában mért QBER alapján, majd minden körben duplázódik
    for p = 1:n_pass
        k_blk(p) = min(n_f, ceil(0.73/Q_k * 2^(p-1))); % a blokkméret nem lehet nagyobb a frame-nél
        perm{p}  = randperm(n_f);       % véletlen sorrend
        pos{p}(perm{p}) = 1:n_f;        % a keverés utáni pozíciók

        n_blk  = ceil(n_f / k_blk(p));  % a blokkok száma, az utolsó blokk rövidebb lehet
        diff_p = [fa(perm{p}) ~= fb(perm{p}), false(1, n_blk*k_blk(p) - n_f)];  % 1, ahol Alice és Bob bitje eltér,  false-al feltöltjük a paritás számoláshoz
        odd{p} = mod(sum(reshape(diff_p, k_blk(p), n_blk), 1), 2) == 1;         % a blokkok paritása (true = páratlan)
        leak   = leak + n_blk;      % minden blokk paritása nyilvános, ezért a leak-be számít

        todo = [repmat(p, 1, nnz(odd{p})); find(odd{p})];    % javítandó (kör, blokk) párok
        % amiíg van javítandó blokk, addig keresünk hibát bináris kereséssel
        while ~isempty(todo)
            % utolsó javítandó blokk: a kör és a blokkindex
            % mindig a legutolsó blokkot javítjuk, mert az új hibák a korábbi blokkokban is lehetnek
            pq  = todo(1,end);
            blk = todo(2,end);
            todo(:,end) = [];
            if ~odd{pq}(blk)
                continue;           % közben páros lett, a következő blokkot javítjuk
            end

            % bináris keresés: a bal fél paritásából a jobb felé is következik
            ib = perm{pq}((blk-1)*k_blk(pq) + 1 : min(blk*k_blk(pq), n_f)); % a javítandó blokk indexei a keverés előtti frame-ben

            % a blokkban páratlan sok hiba van
            while numel(ib) > 1
                fel  = floor(numel(ib) / 2);
                % Alice elküldi a bal fél paritását, Bob kiszámolja a sajátját, és összevetik
                leak = leak + 1;
                % Ha a bal fél páratlan, akkor a hiba ott van, különben a jobb félben. A bal fél indexei: ib(1:fel)
                if mod(nnz(fa(ib(1:fel)) ~= fb(ib(1:fel))), 2) == 1
                    ib = ib(1:fel);
                else
                    ib = ib(fel+1:end);
                end
            end
            % a blokkban egyetlen hiba maradt, azt kijavítjuk
            fb(ib) = ~fb(ib);

            % a javított bit minden eddigi kör blokkjának paritását megfordítja
            for r = 1:p
                % a javított bit a r-edik kör keverése után melyik blokkba esik
                jb = ceil(pos{r}(ib) / k_blk(r));
                % a blokk paritása megfordul, ha eddig páratlan volt, most páros lesz, és fordítva
                odd{r}(jb) = ~odd{r}(jb);
                if odd{r}(jb)
                    % a blokk most páratlan, ezért újra kell keresni a hibát bináris kereséssel
                    todo(:,end+1) = [r; jb]; %#ok<SAGROW>
                end
            end
        end
    end
    % kijavított frame Bobnál
    kb(idx) = fb;

    % --- Megerősítés: mindketten kiszámolják ugyanannak a véletlen hash-nek a képét ---
    % Kétuniverzális Toeplitz-hash: y = T x (mod 2)
    % T egy t_conf x n_f-es véletlen Toeplitz-mátrix, T(i,j) = seed(i - j + n_f)
    
    % A hash hossza t_conf bit, a frame hossza n_f bit
    % A seed hossza n_f + t_conf - 1 bit, a Toeplitz-mátrix minden átlója azonos
    % A toeplitz(c, r) az első oszlopból (c) és az első sorból (r) építi fel a mátrixot:
    %   c = seed(n_f : n_f+t_conf-1)   első oszlop, j = 1:  T(i,1) = seed(i - 1 + n_f)
    %   r = seed(n_f : -1 : 1)         első sor,    i = 1:  T(1,j) = seed(1 - j + n_f)

    % véletlen bitet sorsolunk a seed-hez, a seed hossza n_f + t_conf - 1, Alice és Bob ugyanazt használja, a seed nyilvános
    % Minden frame-hez újat sorsol, a Cascade után
    seed = rand(1, n_f + t_conf - 1) < 0.5;
    % A Toeplitz-mátrixot a seed-ből építjük fel, majd a frame-eket hash-eljük
    T    = toeplitz(double(seed(n_f : n_f+t_conf-1)), double(seed(n_f : -1 : 1)));
    % A hash-elt kulcsok: Alice és Bob frame-jeit hash-eljük, majd összevetjük
    hash = mod(T * double([fa(:), fb(:)]), 2) == 1;   % 1. oszlop: Alice hash-e, 2. oszlop: Bobé
    if isequal(hash(:,1), hash(:,2))
        kept(idx) = true;                   % a frame bitjei maradnak
        leak_kept = leak_kept + leak;       % a nyilvános paritásbitek száma a megtartott frame-ekből
        n_kept    = n_kept + 1;             % a megtartott frame-ek száma
    end
end
% Alice kulcsa, Bob kulcsa, a megerősítésen átment bitek, a nyilvános paritásbitek száma
ka = ka(kept);
kb = kb(kept);
n_rec = numel(ka);

fprintf('2. Cascade: %d frame, %d kör, első blokkméret %d bit | %d hiba volt a maradékban\n', ...
    n_frames, n_pass, min(n_frame, ceil(0.73/Q_k)), n_err_pre);
fprintf('3. Megerősítés (%d bites hash): %d/%d frame maradt, %d bit | észrevétlen hiba: %d\n', ...
    t_conf, n_kept, n_frames, n_rec, nnz(ka ~= kb));
if n_rec == 0
    fprintf('   Egy frame sem ment át a megerősítésen, nincs kulcs.\n');
    return;
end

% ----- Egyeztetett kulcs (titkosságerősítés nélkül) -----
% A hibajavítás és a megerősítés utáni közös kulcs. Ez még NEM titkos: Eve ismerheti
% egy részét (kvantumcsatorna) és látta a nyilvános paritás- és hash-biteket.
T_real = N_total / f_rep;
fprintf('\n--- Egyeztetett kulcs (titkosságerősítés nélkül) ---\n');
fprintf('Siftelt %d -> becslés után %d -> megerősítés után %d bit\n', n_sift, n_pp, n_rec);
fprintf('Egyeztetett kulcsráta: %.4g bit/s | nyilvános bitek: %d paritás + %d hash (%.2f%%)\n', ...
    n_rec/T_real, leak_kept, t_conf*n_kept, 100*(leak_kept + t_conf*n_kept)/n_rec);
fprintf('Alice és Bob kulcsa azonos: %d | Bias: P(0) = %.6f\n', isequal(ka, kb), mean(~ka));
fprintf('Utófeldolgozás futásideje: %.1f s\n', toc);