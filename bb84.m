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

clear; clc; rng('shuffle');

%% Paraméterek
f_rep        = 1e6;     % órajel [Hz]
T_sim        = 1;       % szimulált idő [s]
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
