% --- BBM92 szimuláció: forrás -> zajos, csillapító csatornák -> Alice és Bob (detektorokkal, afterpulse) ---
% Sűrűségmátrixos leírás, ideális összefonódott forrás, mindkét karon bit-flip és
% depolarizáló zaj + csillapítás, Alice-nál és Bobnál is nem ideális, afterpulse-os detektorok.
%
%   Forrás:   minden pulzusban egy összefonódott pár, |Phi+> = (|00> + |11>)/sqrt(2),
%             rho = |Phi+><Phi+| (4x4); az egyik qubit Alice-hoz, a másik Bobhoz megy
%   Csatorna: karonként Kraus-operátorok, ugyanaz a zaj, mint a BB84-ben:
%               depolarizáció: rho -> (1-p) rho + p/3 (X rho X + Y rho Y + Z rho Z)
%               bit-flip:      rho -> (1-p) rho + p X rho X
%             a közös állapotra kron(K_A, K_B) hat (a csillapítás a POVM-ben van)
%   Mérés:    POVM oldalanként, a csatorna csillapításával együtt:
%             E_k = T * eta_k * Pi_k  (k = D0, D1),  E_nincs = I - E_0 - E_1
%             P(x, y) = Tr((E_x (x) E_y) rho)   ((x) = tenzorszorzat, MATLAB-ban kron)
%             x, y = D0, D1 vagy nincs jelklikk -> 3 x 3 = 9 kimenet bázispáronként
%   Detektorok, mindkét oldalon két detektor:
%             jitter, sötétzaj, kapuzás, holtidő, dupla klikk -> eldobva
%   Afterpulse, mind a 4 detektoron:
%             egy lavina p_ap eséllyel csapdába ejt egy töltéshordozót, amely
%             Exp(tau_ap) idő múlva szabadul fel. Ha ez a holtidő után, egy kapun
%             belül történik, ugyanaz a detektor hamis klikket ad, ami újabb
%             afterpulse-t indíthat.
%   Sifting:  csak a koincidenciás (mindkét oldalon pontosan egy klikk) ÉS egyező
%             bázisú párok maradnak
%   Utófeldolgozás

clear; clc; rng('shuffle');

%% Paraméterek
f_rep        = 1e6;     % órajel [Hz] (1 MHz-en a következő kapu 1 us múlva nyílik,
                        % addigra a tau_ap = 100 ns-os afterpulse lecseng)
T_sim        = 10;       % szimulált idő [s]
chunk_pulses = 1e6;     % egy blokk mérete (memória/sebesség tradeoff)

dt       = 1 / f_rep;   % pulzusperiódus [s]
N_total  = round(T_sim * f_rep);
N_chunks = ceil(N_total / chunk_pulses);

%% Csatorna — paraméter
L_A   = 10;             % forrás -> Alice [km]
L_B   = 10;             % forrás -> Bob [km]
alpha = 0.2;            % [dB/km]
T_A   = 10^(-(alpha*L_A)/10); % Alce felé csillapítás: T = 10^(-alpha*L/10) = 10^(-dB/10)
T_B   = 10^(-(alpha*L_B)/10); % Bob felé csillapítás: T = 10^(-alpha*L/10) = 10^(-dB/10)

p_bf  = [0.02, 0.02];   % bit-flip valószínűsége: forrás -> Alice, forrás -> Bob
p_dep = [0.03, 0.03];   % depolarizáció valószínűsége: forrás -> Alice, forrás -> Bob

%% Detektorok — paraméter (sor: 1 = Alice, 2 = Bob; oszlop: D0 = bit 0, D1 = bit 1)
eta_D       = [0.10, 0.10;      % detektor-hatásfok: Alice
               0.10, 0.10];     %                    Bob
lambda_dark = [200, 200;        % sötétszámlálási ráta [Hz]: Alice
               200, 200];       %                            Bob
tau_dead    = 50e-9;            % holtidő [s] (nem bénuló)
T_gate      = 200e-12;          % detekciós kapu szélessége [s]
sigma_jit   = 30e-12;           % időbeli jitter (1 sigma) [s]
p_ap        = [0.05, 0.05;      % P(egy lavina csapdába ejt egy töltéshordozót): Alice
               0.05, 0.05];     %                                                Bob
tau_ap      = 100e-9;           % a csapdázott töltéshordozó élettartama [s]

assert(T_gate <= dt, 'A kapu (T_gate) nem lehet hosszabb a pulzusperiódusnál.');

%% Utófeldolgozás — paraméter
v_pe    = 0.5;          % a paraméterbecslésre feláldozott hányad (előadás: BB84-nél a fele)
Q_tol   = 0.11;         % QBER-küszöb: e fölött a protokoll leáll
eps_pe  = 1e-10;        % P(a maradék bitek QBER-je meghaladja a becsült felső korlátot)
n_frame = 1e4;          % a Cascade egy frame-jének hossza [bit]
n_pass  = 4;            % Cascade körök száma
t_conf  = 64;           % a megerősítő hash hossza [bit]

%% Állapot (forrás)
ket0 = [1; 0];
ket1 = [0; 1];
ketP = (ket0 + ket1) / sqrt(2);
ketM = (ket0 - ket1) / sqrt(2);

% Összefonódott pár: a két qubit közös állapota a kron-nal áll össze (4 dimenziós tér:
% |00>, |01>, |10>, |11>; az első qubit Alice-é, a második Bobé)
psi = (kron(ket0, ket0) + kron(ket1, ket1)) / sqrt(2);    % |Phi+>
rho = psi * psi';                                         % 4x4 sűrűségmátrix

% psi = [1; 0; 0; 1] / sqrt(2)
% rho = 1/2 · [1 0 0 1
%              0 0 0 0
%              0 0 0 0
%              1 0 0 1]

%% Mérés (Alice és Bob)
% Ugyanazok a projektorok, mint a BB84-ben, és mindkét oldal ugyanezeket használja:
% Pi(:,:,kimenet,bázis); kimenet 1 = D0 (bit 0), 2 = D1 (bit 1)
%                        bázis   1 = Z,          2 = X
%                        Z: |0><0|, |1><1|; X: |+><+|, |-><-|
Pi = zeros(2, 2, 2, 2);
Pi(:,:,1,1) = ket0 * ket0';  Pi(:,:,2,1) = ket1 * ket1';
Pi(:,:,1,2) = ketP * ketP';  Pi(:,:,2,2) = ketM * ketM';

% POVM oldalanként: E(:,:,kimenet,bázis,oldal); kimenet 3 = nincs jelklikk
%   E_k = T * eta_k * Pi_k,  E_nincs = I - E_0 - E_1
% A csillapítás azért került a mérésbe a Kraus-operátor helyett, mert most az is
% számít, ha csak az EGYIK foton érkezik meg: a másik oldalon egy sötétklikk ilyenkor
% hamis koincidenciát adhat. A kron(sqrt(T_A) I, sqrt(T_B) I) Kraus-operátor csak azt
% az esetet írná le, amikor mindkét foton megérkezik. Mivel a veszteség nem függ a
% polarizációtól, "elveszett a szálban" és "megérkezett, de a detektor nem érzékelte"
% ugyanazt jelenti: nincs jelklikk. Így a 9 kimenet valószínűsége összesen 1
% Az E tömb 5 dimenziós: E(:,:,kimenet,bázis,oldal), ahol
% kimenet 1 = D0, 2 = D1, 3 = nincs jelklikk
% bázis 1 = Z, 2 = X
% oldal 1 = Alice, 2 = Bob
T_side = [T_A, T_B];
E = zeros(2, 2, 3, 2, 2);
for s = 1:2
    for b = 1:2
        E(:,:,1,b,s) = T_side(s) * eta_D(s,1) * Pi(:,:,1,b);
        E(:,:,2,b,s) = T_side(s) * eta_D(s,2) * Pi(:,:,2,b);
        E(:,:,3,b,s) = eye(2) - E(:,:,1,b,s) - E(:,:,2,b,s);
    end
end

%% Csatorna (Kraus-operátorok)
I2 = eye(2);
X  = [0 1; 1 0];                    % Pauli-mátrixok
Y  = [0 -1i; 1i 0];
Z  = [1 0; 0 -1];

% Oldalanként (1 = forrás -> Alice, 2 = forrás -> Bob) ugyanaz, mint a BB84-ben:
% depolarizáció, majd bit-flip; egymás után alkalmazva a Kraus-operátorok minden
% párosításának szorzata, Oldalanként 4 x 2 = 8 operátor.
K_kar = cell(1, 2);
for s = 1:2
    K_dep = {sqrt(1 - p_dep(s)) * I2, sqrt(p_dep(s)/3) * X, sqrt(p_dep(s)/3) * Y, sqrt(p_dep(s)/3) * Z};
    K_bf  = {sqrt(1 - p_bf(s)) * I2, sqrt(p_bf(s)) * X};
    K_kar{s} = {};
    for a2 = 1:numel(K_bf)
        for a3 = 1:numel(K_dep)
            K_kar{s}{end+1} = K_bf{a2} * K_dep{a3};
        end
    end
end

% A két oldal egymástól független, ezért a közös 4x4-es állapotra ható Kraus-operátorok
% a két oldal operátorainak tenzorszorzatai: kron(K_A{i}, K_B{j}), itt 8 x 8 = 64 darab.
K = {};
for i = 1:numel(K_kar{1})
    for j = 1:numel(K_kar{2})
        K{end+1} = kron(K_kar{1}{i}, K_kar{2}{j}); %#ok<SAGROW>
    end
end

% Ellenőrzés: sum_k K_k' K_k = I (a zaj nyomtartó; a csillapítás a POVM-ben van)
KK = zeros(4);
for k = 1:numel(K)
    KK = KK + K{k}' * K{k};
end

rho_B = zeros(4, 4);
for k = 1:numel(K)
    rho_B = rho_B + K{k} * rho * K{k}';
end

%% Táblázat: P(a, b, x, y) = Tr((E_x^a,Alice (x) E_y^b,Bob) * rho_B)
% Bázispárenként (sor: index = a + 2*(b-1)) a 9 kimenet (oszlop: o = x + 3*(y-1)),
% ahol x, y = 1: D0, 2: D1, 3: nincs jelklikk. Ez a jelkliikkek eloszlása; a
% sötétzaj és a holtidő a detektor-részben jön
% E(:,:,x,a,1): Alice (5. index = 1) POVM-eleme, x kimenet, a bázis
% E(:,:,y,b,2): Bob (5. index = 2) POVM-eleme, y kimenet, b bázis
Pt = zeros(4, 9);
for a = 1:2
    for b = 1:2
        for x = 1:3
            for y = 1:3
                Pt(a + 2*(b-1), x + 3*(y-1)) = ...
                    real(trace(kron(E(:,:,x,a,1), E(:,:,y,b,2)) * rho_B));   % számitásnál kerekítési hiba miatt kis képzetes rész lehet, ezért real() függvény
            end
        end
    end
end
C = cumsum(Pt, 2);                  % kumulált valószínűségek a sorsoláshoz

% |─o1-|o2|----o3-----|o4|-o5-|----o6-----|----o7----|----o8---|---------o9----------|

fprintf('Csatorna: L_A = %g km, L_B = %g km, alpha = %g dB/km -> T_A = %.3g, T_B = %.3g\n', ...
    L_A, L_B, alpha, T_A, T_B);
fprintf('Zaj (Alice | Bob oldal): bit-flip p = [%g | %g], depolarizáció p = [%g | %g]\n', ...
    p_bf(1), p_bf(2), p_dep(1), p_dep(2));
fprintf('  %d Kraus-operátor, sum K''K = I eltérése: %.1g\n', numel(K), max(abs(KK - eye(4)), [], 'all'));
fprintf('Detektorok (Alice | Bob): eta = [%g %g | %g %g], sötét = [%g %g | %g %g] Hz\n', ...
    eta_D(1,:), eta_D(2,:), lambda_dark(1,:), lambda_dark(2,:));
fprintf('  holtidő = %g ns, kapu = %g ps, jitter = %g ps\n', tau_dead*1e9, T_gate*1e12, sigma_jit*1e12);
fprintf('Afterpulse (Alice | Bob): p_ap = [%g %g | %g %g], tau_ap = %g ns\n', ...
    p_ap(1,:), p_ap(2,:), tau_ap*1e9);
bazis = {'Z', 'X'};
fprintf('\nJelklikk-valószínűségek páronként (A = csak Alice-nál, B = csak Bobnál):\n');
fprintf('Bázis (A,B) | P(00)      P(01)      P(10)      P(11)      | P(A)       P(B)       P(egyik sem)\n');
for a = 1:2
    for b = 1:2
        r = Pt(a + 2*(b-1), :);
        fprintf('    %s, %s     | %.3g  %.3g  %.3g  %.3g  | %.3g  %.3g  %.3g\n', bazis{a}, bazis{b}, ...
            r(1), r(4), r(2), r(5), r(3) + r(6), r(7) + r(8), r(9));
    end
end

%% Akkumulátorok
% A sifting és QBER számításához
% tot_sifted: összes siftelt bit
% tot_errors: összes siftelt bitből hány hibás
% tot_n0: összes siftelt bitből hány 0
tot_detected = 0; tot_sifted = 0; tot_errors = 0; tot_n0 = 0;
tot_sifted_b = [0 0];               % siftelt bitek bázisonként (Z, X)
tot_errors_b = [0 0];               % hibák bázisonként (Z, X)
tot_single   = [0 0];               % oldalanként az egy klikkes pulzusok (Alice, Bob)
tot_double   = [0 0];               % oldalanként a dupla klikkek (eldobva)
tot_src      = zeros(4, 3);         % detektor x forrás (1 = jel, 2 = sötét, 3 = afterpulse)
last_click   = -inf(1, 4);          % holtidő-állapot: utolsó klikk abszolút ideje
ap_pending   = repmat({zeros(1,0)}, 1, 4);   % függő felszabadulások abszolút ideje detektoronként
det_nev      = {'Alice D0', 'Alice D1', 'Bob D0', 'Bob D1'};
key_a        = cell(1, N_chunks);   % Alice és Bob siftelt kulcsa chunkonként
key_b        = cell(1, N_chunks);

fprintf('\nIndulás: %d pár (%g s @ %g GHz), %d chunk\n', ...
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

    % ----- Bázisválasztás (egymástól függetlenül) -----
    alice_bases = randi([0 1], 1, N);           % 0 = Z, 1 = X
    bob_bases   = randi([0 1], 1, N);
    index = alice_bases + 1 + 2*bob_bases;      % a táblázat sora
    % index = 1 -> (Z,Z), index = 2 -> (X,Z), index = 3 -> (Z,X), index = 4 -> (X,X) Alice, Bob bázisok sorrendjében


    % ----- Jelklikkek: a két oldal együtt, egy véletlen számmal, 9 kimenet -----
    % u a kumulált sávok közül abba esik, amelyik előtt C(index, o) <= u;
    %  a sorszáma o = 1 + (hány határ van u alatt). Ebből x = Alice, y = Bob kimenete.
    u = rand(N, 1);

    % C(index, 1:8) = minden pulzushoz kiveszi a saját bázispárjának 8 határát
    % U>=C(index, 1:8) logikai mátrix, soronként a 8 határhoz tartozó igaz/hamis értékek
    % sum(..., 2) = soronként összegzi, hány határ van u alatt, ebből szakasz sorszámot csinál
    o = 1 + sum(u >= C(index, 1:8), 2)';
    
    % o − 1 = (x−1) + 3·(y−1)
    x = mod(o - 1, 3) + 1;                      % 1 = D0, 2 = D1, 3 = nincs jelklikk, Alice
    y = floor((o - 1) / 3) + 1;                 % 1 = D0, 2 = D1, 3 = nincs jelklikk, Bob
    sig_idx = {find(x == 1), find(x == 2), find(y == 1), find(y == 2)};   % szétválogatás detektoronként

    % ----- Detektorok: jitter, sötétzaj, kapuzás, holtidő, afterpulse (4 detektor) -----
    pulse_clicks = false(4, N);                 % sor: Alice D0, Alice D1, Bob D0, Bob D1
    for d = 1:4
        s  = ceil(d / 2);                       % oldal: 1 = Alice, 2 = Bob
        dd = d - 2*(s - 1);                     % detektor az oldalon: 1 = D0, 2 = D1

        % jel: a k-adik pulzus a periódus közepén, (k - 1/2)*dt-kor érkezik, + jitter
        t_sig  = (sig_idx{d} - 0.5) * dt + sigma_jit * randn(1, numel(sig_idx{d}));
        % sötétzaj: Poisson-folyamat a teljes chunkra (a kapun kívül is)
        % a darabszám Poisson-eloszlású, az időpontok egyenletesen szórnak
        % itt kell a Statistics and Machine Learning Toolbox, a poissrnd() függvény miatt
        N_dark = poissrnd(lambda_dark(s, dd) * T_chunk);
        t_dark = T_chunk * rand(1, N_dark);

        % a jel és sötétzaj időpontjait egyesítjük, és sorba rendezzük
        % a forrás típusát is tároljuk (1 = jel, 2 = sötét)
        t   = [t_sig, t_dark];
        src = [ones(1, numel(t_sig)), 2*ones(1, numel(t_dark))];

        % kapuzás: csak a pulzusidő körüli T_gate ablakba eső esemény számít
        k = round(t/dt + 0.5);                  % a legközelebbi pulzus indexe
        in_gate = k >= 1 & k <= N & abs(t - (k - 0.5)*dt) <= T_gate/2;
        t = t(in_gate);  src = src(in_gate);  k = k(in_gate);
        % sorba rendezzük az eseményeket idő szerint, hogy a holtidőt alkalmazhassuk
        [t, ord] = sort(t);
        src = src(ord);  k = k(ord);

        % holtidő + afterpulse: egy klikk után tau_dead ideig a detektor vak.
        % előző klikk ideje: last_click(d), a chunkon belüli idő: last_click(d) - t0
        % Időrendben kell haladni, mert minden elfogadott klikk új afterpulse-t
        % indíthat, és az afterpulse is klikk: holtidőt és újabb afterpulse-t okoz.
        % A holtidő-állapot és a függő felszabadulások a chunkok közt is megmaradnak.
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
                if rand < p_ap(s, dd)
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

    % ----- Bitek oldalanként -----
    %   csak D0 -> 0 | csak D1 -> 1 | egyik sem -> -1
    %   mindkettő (dupla klikk) -> a bitet eldobjuk (-1 marad), de számon tartjuk
    bits = -ones(2, N);                         % sor: 1 = Alice, 2 = Bob
    for s = 1:2
        c0 = pulse_clicks(2*s - 1, :);
        c1 = pulse_clicks(2*s, :);
        bits(s,  c0 & ~c1) = 0;
        bits(s, ~c0 &  c1) = 1;
        tot_double(s) = tot_double(s) + nnz(c0 & c1);
        tot_single(s) = tot_single(s) + nnz(bits(s,:) >= 0);
    end
    alice_bits = bits(1,:);
    bob_bits   = bits(2,:);
    detected   = (alice_bits >= 0) & (bob_bits >= 0);   % koincidencia

    % ----- Sifting -----
    % Alice és Bob nyilvánosan egyezteti, hol volt klikk és milyen bázisban mértek;
    % csak a koincidenciás ÉS egyező bázisú párok maradnak.
    sift = detected & (alice_bases == bob_bases);
    tot_detected = tot_detected + nnz(detected);
    tot_sifted   = tot_sifted + nnz(sift);
    tot_errors   = tot_errors + nnz(alice_bits(sift) ~= bob_bits(sift));
    tot_n0       = tot_n0 + nnz(bob_bits(sift) == 0);
    key_a{c}     = alice_bits(sift) == 1;
    key_b{c}     = bob_bits(sift) == 1;

    % bázisonként is: a bit-flip csak a Z bázisban okoz hibát
    for b = 1:2
        sb = sift & (alice_bases == b-1);
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

% Várt értékek holtidő nélkül (közelítés: átlagos hatásfok és sötétzaj oldalanként)
% Egy oldalon pontosan egy klikk kétféleképpen lehet:
%   jel:   a foton az egyik detektoron kattan, és a MÁSIKON nincs sötétklikk: p (1-q)
%   sötét: nincs jelklikk, és pontosan egy detektoron van sötétklikk:   (1-p) 2q (1-q)
% (ha a jel mellé a másik detektoron sötétklikk jön, abból dupla klikk lesz -> eldobva)
p_A    = T_A * mean(eta_D(1,:));                % jelklikk Alice-nál / pulzus
p_B    = T_B * mean(eta_D(2,:));                % jelklikk Bobnál / pulzus
q_A    = mean(lambda_dark(1,:)) * T_gate;       % sötétklikk egy detektoron / pulzus
q_B    = mean(lambda_dark(2,:)) * T_gate;
s_A    = p_A * (1 - q_A);                       % egy klikk, a fotontól
s_B    = p_B * (1 - q_B);
n_A    = (1 - p_A) * 2*q_A * (1 - q_A);         % egy klikk, sötétzajtól
n_B    = (1 - p_B) * 2*q_B * (1 - q_B);
c_true = s_A * s_B;                             % valódi koincidencia: mindkét oldalon a pár fotonja
c_acc  = s_A*n_B + n_A*s_B + n_A*n_B;           % véletlen koincidencia: legalább az egyik sötét

% A csatorna-zaj hibája egyező bázisnál, karonként (oszlop: Alice, Bob karja):
q_dep = 2*p_dep/3;                              % depolarizáció: X és Y fordít Z-ben, Y és Z fordít X-ben
e_kar = [p_bf.*(1-q_dep) + (1-p_bf).*q_dep;     % Z: a bit-flip VAGY a depolarizáció fordít
         q_dep];                                % X: a bit-flip nem hat
% A pár bitjei akkor térnek el, ha pontosan az egyik kar fordított (ha mindkettő, újra egyeznek):
e_b = (e_kar(:,1).*(1-e_kar(:,2)) + (1-e_kar(:,1)).*e_kar(:,2))';
R_tau  = f_rep * ([p_A, p_B]/2 + [q_A, q_B]) * tau_dead;

% Afterpulse / klikk telítés nélkül: a felszabadulás a holtidő után, egy kapuba esik
% (a klikket a kapu közepére véve). A láncokkal együtt egy elsődleges klikkre
% p/(1-p) afterpulse jut. Telítésnél a detektor gyakran vak, ezért a mért érték kisebb.
n_g   = floor(tau_dead/dt) : ceil((tau_dead + 30*tau_ap)/dt);    % a klikk utáni kapuk
lo    = max(n_g*dt - T_gate/2, tau_dead);
hi    = max(n_g*dt + T_gate/2, tau_dead);
p_eff = mean(p_ap(:)) * sum(exp(-lo/tau_ap) - exp(-hi/tau_ap));
a_v   = p_eff / (1 - p_eff);                                     % várt
a_m   = sum(tot_src(:,3)) / max(sum(tot_src(:,1:2), 'all'), 1);  % mért

% QBER a MÉRT klikkekből (a holtidő és az afterpulse így már benne van):
% oldalanként a fotontól és a zajtól (sötét + afterpulse) származó klikkek aránya
% pulzusonként. A két oldal fotonjai egymástól függetlenül vesznek el, ezért a valódi
% koincidencia a két oldal fotonklikkjének szorzata; minden más véletlen koincidencia.
ph_m   = [sum(tot_src(1:2,1)), sum(tot_src(3:4,1))] / N_total;
nz_m   = [sum(tot_src(1:2,2:3), 'all'), sum(tot_src(3:4,2:3), 'all')] / N_total;
c_tr_m = ph_m(1) * ph_m(2);
c_ac_m = ph_m(1)*nz_m(2) + nz_m(1)*ph_m(2) + nz_m(1)*nz_m(2);
% A valódi koincidenciák e_b, a véletlenek 1/2 eséllyel hibásak:
QBER_v_b = (e_b*c_tr_m + 0.5*c_ac_m) / (c_tr_m + c_ac_m);
QBER_v   = mean(QBER_v_b);                      % a két bázis egyforma valószínű

fprintf('\n--- BBM92 (zajos csatorna + detektorok + afterpulse) eredmények ---\n');
fprintf('Szimulált idő: %g s | futásidő: %.1f s\n', N_total/f_rep, toc);
fprintf('Elküldött párok: %d\n', N_total);
for d = 1:4
    fprintf('  %s klikkek: %d  [jel = %d, sötét = %d, afterpulse = %d]\n', det_nev{d}, ...
        sum(tot_src(d,:)), tot_src(d,1), tot_src(d,2), tot_src(d,3));
end
fprintf('Afterpulse / elsődleges klikk: %.4f%%  (várt telítés nélkül: %.4f%%)\n', 100*a_m, 100*a_v);
fprintf('Egy klikkes pulzusok: Alice %d, Bob %d  (várt holtidő nélkül: %.0f, %.0f)\n', ...
    tot_single, N_total*(s_A + n_A), N_total*(s_B + n_B));
fprintf('Dupla klikkek (eldobva): Alice %d, Bob %d\n', tot_double);
fprintf('Koincidenciák: %d  (várt holtidő nélkül: %.0f | a mért klikkekből: %.0f, ebből véletlen: %.0f)\n', ...
    tot_detected, N_total*(c_true + c_acc), N_total*(c_tr_m + c_ac_m), N_total*c_ac_m);
fprintf('  R * tau_dead = %.3g (Alice), %.3g (Bob) detektoronként (ha nem << 1, a holtidő sok klikket elvesz)\n', R_tau);
fprintf('Siftelt kulcs: %d bit | siftelt kulcsráta: %.4g bit/s\n', ...
    tot_sifted, tot_sifted/(N_total/f_rep));
fprintf('QBER:   %.4f%%  (várt a mért klikkekből: %.4f%%)\n', 100*QBER, 100*QBER_v);
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