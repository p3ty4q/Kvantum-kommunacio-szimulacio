function fig = bb84_gui()
% BB84_GUI  Grafikus felület a bb84.m szimulációhoz
%
%   bb84_gui            megnyitja a felületet
%   fig = bb84_gui;     ugyanez, és visszaadja az ablakot (uifigure)
%
% Bal oldalt a paraméterek
% A "Szimuláció indítása" gomb a megadott paraméterekkel lefuttatja a bb84.m-et;
% jobb oldalt megjelenik az eredménytáblázat és a szkript kiírása.
%
% A felület az alapértékeket a bb84.m-ből olvassa ki: minden paraméternél az
% első "név = érték;" sort keresi meg.
% Indításkor a bb84.m szövegében ugyanezekben a sorokban kicseréli az értéket a
% mezőben megadottra, és ezt a másolatot futtatja le. A másolat
% ideiglenes fájl a bb84.m mellett; a futás után törlődik.
% A bb84.m-hez így nem kell semmit hozzáírni: a clear, a származtatott
% mennyiségek (dt, N_total, transmittance, ...) ugyanolyan, mint önálló futtatáskor.
%
% Az ablak programból készül (uifigure + uigridlayout).

P         = parameter_tabla();
bb84_fajl = fullfile(fileparts(mfilename('fullpath')), 'bb84.m');
alap      = alapertekek(bb84_fajl, P(:,1));

fig = uifigure('Name', 'BB84 szimuláció', 'Position', [60 60 1350 820]);
fo  = uigridlayout(fig, [1 2], 'ColumnWidth', {520, '1x'});   % bal oldal: felirat (név) + mezők

%% Bal oldal: paraméterek
bal   = uigridlayout(fo, [5 1], 'RowHeight', repmat({'fit'}, 1, 5), ...
    'Scrollable', 'on', 'Padding', [0 0 4 0]);
cimek = {'Forrás és futás', 'Csatorna', 'Bob detektorai', 'Utófeldolgozás'};
mezo  = struct();                       % mezo.(név): a paraméter beviteli mezője(i)
for pn = 1:numel(cimek)
    sorok  = find([P{:,7}] == pn);
    vektor = any(arrayfun(@(i) numel(alap.(P{i,1})) > 1, sorok));   % van-e D0/D1 oszlop
    n_sor  = numel(sorok) + vektor;
    g = uigridlayout(uipanel(bal, 'Title', cimek{pn}, 'FontWeight', 'bold'), [n_sor 4], ...
        'ColumnWidth', {'1x', 85, 85, 45}, 'RowHeight', repmat({22}, 1, n_sor), ...
        'RowSpacing', 3, 'Padding', [8 5 8 5]);
    sor = 1;
    if vektor                           % fejléc a detektoronkénti paraméterekhez
        cimke(g, 'D0 (bit 0)', sor, 2, 'center');
        cimke(g, 'D1 (bit 1)', sor, 3, 'center');
        sor = sor + 1;
    end
    for i = sorok
        nev = P{i,1};
        ert = alap.(nev) / P{i,4};      % a felületen mutatott mértékegységben
        cimke(g, sprintf('%s (%s)', P{i,2}, nev), sor, 1, 'left');   % felirat (név a bb84.m-ben)
        if P{i,6}, formatum = '%.0f'; else, formatum = '%.4g'; end   % egész: 10000, nem 1e+04
        h = matlab.ui.control.NumericEditField.empty;
        for d = 1:numel(ert)            % detektoronkénti paraméternél két mező
            h(d) = uieditfield(g, 'numeric', 'Value', kerekit(ert(d)), 'Limits', P{i,5}, ...
                'RoundFractionalValues', P{i,6}, 'ValueDisplayFormat', formatum, ...
                'ValueChangedFcn', @(~, ~) varhato_frissit(fig));
            h(d).Layout.Row    = sor;
            h(d).Layout.Column = 1 + d;
        end
        if isscalar(ert), h.Layout.Column = [2 3]; end
        cimke(g, P{i,3}, sor, 4, 'left');
        mezo.(nev) = h;
        sor = sor + 1;
    end
end

% Várható értékek: a paraméterekből, szimuláció nélkül (élőben frissül)
gi   = uigridlayout(uipanel(bal, 'Title', 'Várható értékek', ...
    'FontWeight', 'bold'), [2 1], 'RowHeight', {'fit', 'fit'}, 'RowSpacing', 2, 'Padding', [8 6 8 6]);
info = uilabel(gi, 'Text', '', 'VerticalAlignment', 'top');
figy = uilabel(gi, 'Text', '', 'WordWrap', 'on', 'FontColor', [208 59 59]/255);

%% Jobb oldal: gombok és állapot, eredménytáblázat, a szkript kiírása
jobb    = uigridlayout(fo, [3 1], 'RowHeight', {44, 445, '1x'}, 'Padding', [4 0 0 0]);
fent    = uigridlayout(jobb, [1 3], 'ColumnWidth', {170, 190, '1x'}, 'RowHeight', {'1x'}, ...
    'Padding', [0 0 0 0]);
gomb    = uibutton(fent, 'Text', 'Szimuláció indítása', 'FontWeight', 'bold', ...
    'BackgroundColor', [42 120 214]/255, 'FontColor', [1 1 1], ...
    'ButtonPushedFcn', @(~, ~) inditas(fig));
uibutton(fent, 'Text', 'Alapértékek a bb84.m-ből', 'ButtonPushedFcn', @(~, ~) alap_visszaallit(fig));
allapot = uilabel(fent, 'Text', 'Add meg a paramétereket, majd indítsd el a szimulációt.', ...
    'FontSize', 13, 'FontWeight', 'bold', 'WordWrap', 'on');
tabla   = uitable(jobb, 'ColumnName', {'Mennyiség'; 'Érték'}, 'RowName', {}, ...
    'ColumnWidth', {'1x', 200}, 'Data', cell(0, 2));
gn    = uigridlayout(uipanel(jobb, 'Title', 'A bb84.m kiírása', 'FontWeight', 'bold'), [1 1], ...
    'Padding', [4 4 4 4]);
naplo = uitextarea(gn, 'Editable', 'off', 'FontName', 'Consolas', ...
    'Value', {'A futás után itt jelenik meg a szkript teljes konzolkiírása.'});

% Az ablak állapota: a visszahívások (callback) innen érik el a vezérlőket.
u = struct('P', {P}, 'mezo', mezo, 'bb84_fajl', bb84_fajl, 'info', info, 'figy', figy, ...
    'allapot', allapot, 'tabla', tabla, 'naplo', naplo, 'gomb', gomb);
fig.UserData = u;
varhato_frissit(fig);

if nargout == 0, clear fig; end
end


function P = parameter_tabla()
% A felületen megjelenő paraméterek, oszloponként:
%   1 név a bb84.m-ben      2 felirat           3 mértékegység a felületen
%   4 skála (bb84.m-beli érték = kijelzett érték * skála)
%   5 [min max] a kijelzett értékre   6 egész szám-e   7 panel (1..4)
% Hogy egy paraméter detektoronkénti-e (D0, D1), azt a bb84.m-beli alapértéke dönti el.
P = {
    'f_rep',       'Órajel',                       'MHz',    1e6,   [0 Inf],  false, 1
    'T_sim',       'Szimulált idő',                's',      1,     [0 Inf],  false, 1
    'L',           'Távolság',                     'km',     1,     [0 Inf],  false, 2
    'alpha',       'Csillapítás',                  'dB/km',  1,     [0 Inf],  false, 2
    'p_bf',        'Bit-flip valószínűsége',       '',       1,     [0 1],    false, 2
    'p_dep',       'Depolarizáció valószínűsége',  '',       1,     [0 1],    false, 2
    'eta_D',       'Detektor-hatásfok',            '',       1,     [0 1],    false, 3
    'lambda_dark', 'Sötétszámlálási ráta',         'Hz',     1,     [0 Inf],  false, 3
    'p_ap',        'Afterpulse valószínűsége',     '',       1,     [0 1],    false, 3
    'tau_dead',    'Holtidő',                      'ns',     1e-9,  [0 Inf],  false, 3
    'T_gate',      'Kapuszélesség',                'ps',     1e-12, [0 Inf],  false, 3
    'sigma_jit',   'Jitter (1 szigma)',            'ps',     1e-12, [0 Inf],  false, 3
    'tau_ap',      'Csapda élettartama',           'ns',     1e-9,  [0 Inf],  false, 3
    'v_pe',        'Mintavételi arány',            '',       1,     [0 1],    false, 4
    'Q_tol',       'QBER-küszöb',                  '%',      1e-2,  [0 50],   false, 4
    'eps_pe',      'Becslés hibavalószínűsége',    '',       1,     [0 1],    false, 4
    'n_frame',     'Frame-hossz',                  'bit',    1,     [1 Inf],  true,  4
    'n_pass',      'Cascade-körök',                '',       1,     [1 Inf],  true,  4
    't_conf',      'Hash-hossz',                   'bit',    1,     [1 Inf],  true,  4
    };
end


function alap = alapertekek(bb84_fajl, nevek)
% Az alapértékek a bb84.m-ből: az első "név = érték;" kezdetű sor (megjegyzés nélkül)
szoveg = fileread(bb84_fajl, 'Encoding', 'UTF-8');
alap   = struct();
for i = 1:numel(nevek)
    tok  = regexp(szoveg, param_minta(nevek{i}), 'tokens', 'once', 'lineanchors');
    ertek = [];
    if ~isempty(tok)
        ertek = str2num(tok{2}); %#ok<ST2NM> % a [0.10, 0.10] alakú vektorok miatt kell
    end
    if isempty(ertek)
        error('bb84_gui:alapertek', 'Nem találom a(z) %s alapértékét a %s fájlban.', ...
            nevek{i}, bb84_fajl);
    end
    alap.(nevek{i}) = ertek;
end
end


function minta = param_minta(nev)
% Egy paraméter sora a bb84.m-ben: "név = érték;" (a sor elején, utána jöhet megjegyzés)
%   1. csoport: "név = "    2. csoport: az érték    3. csoport: ";"
minta = ['^([ \t]*' nev '[ \t]*=[ \t]*)([^;%\n]+)(;)'];
end


function par = parameterek(u)
% A mezők értékei SI-egységekben, a bb84.m változóneveivel: ez lesz a gui_par
par = struct();
for i = 1:size(u.P, 1)
    par.(u.P{i,1}) = [u.mezo.(u.P{i,1}).Value] * u.P{i,4};
end
end


function varhato_frissit(fig)
% A paraméterekből számolt várható értékek, a bb84.m "várt" képleteivel
u      = fig.UserData;
par    = parameterek(u);
T      = 10^(-par.alpha * par.L / 10);                                   % átvitel
N      = round(par.T_sim * par.f_rep);                                   % pulzusok
p_prim = T * mean(par.eta_D) + sum(par.lambda_dark) * par.T_gate;       % jel + sötétklikk / pulzus
n_sift = N * p_prim / 2;                % várható siftelt bitek, csak a figyelmeztetéshez
u.info.Text = {
    sprintf('Csatorna: %.3g dB veszteség, átvitel T = %.3g', par.alpha * par.L, T)
    sprintf('Pulzusok: %s  (%g s × %g MHz)', ezres(N), par.T_sim, par.f_rep / 1e6)
    };
if par.f_rep > 0 && par.T_gate > 1 / par.f_rep
    u.figy.Text = 'Hiba: a kapu (T_gate) hosszabb, mint a pulzusperiódus (1/f_rep).';
elseif n_sift < 1e4
    u.figy.Text = 'Figyelem: kevés siftelt bit (< 10 000), a QBER-becslés pontatlan lesz. Növeld a szimulált időt.';
else
    u.figy.Text = '';
end
end


function inditas(fig)
% "Szimuláció indítása": a bb84.m futtatása a mezők értékeivel
varhato_frissit(fig);
u   = fig.UserData;
par = parameterek(u);
if par.f_rep <= 0 || par.T_sim <= 0
    uialert(fig, 'Az órajel és a szimulált idő legyen pozitív.', 'Hibás paraméter');
    return;
end
if par.T_gate > 1 / par.f_rep
    uialert(fig, 'A kapu (T_gate) nem lehet hosszabb a pulzusperiódusnál (1/f_rep).', 'Hibás paraméter');
    return;
end

N = round(par.T_sim * par.f_rep);
u.gomb.Enable          = 'off';
u.allapot.Text         = 'A szimuláció fut...';
u.allapot.FontColorMode = 'auto';      % a téma (világos/sötét) szerinti alapszín
dlg = uiprogressdlg(fig, 'Title', 'BB84 szimuláció', 'Indeterminate', 'on', ...
    'Message', sprintf('%s pulzus szimulálása...', ezres(N)));
drawnow;

t0 = tic;
try
    [ered, naplo] = bb84_futtat(par, u.bb84_fajl);
catch hiba
    close(dlg);
    u.gomb.Enable       = 'on';
    u.allapot.Text      = 'Hiba: a bb84.m futása megszakadt.';
    u.allapot.FontColor = [208 59 59]/255;
    % az ideiglenes másolat neve helyett a bb84.m-et írja ki (a sorszámok azonosak)
    uialert(fig, regexprep(hiba.message, 'bb84_gui_tp\w*', 'bb84'), 'Hiba a bb84.m futása közben');
    return;
end
ido = toc(t0);
close(dlg);
u.gomb.Enable = 'on';
eredmenyek(u, ered, naplo, ido);
end


function [ered, naplo] = bb84_futtat(par, bb84_fajl)
% A bb84.m futtatása a felület paramétereivel.
% A bb84.m szövegében minden paraméter első "név = érték;" sorában az értéket a
% felületére cseréli (pl. "L = 20;" -> "L = 35;"), a sorok száma nem változik.
% Ezt ideiglenes fájlba írja a bb84.m mellé (így a mappára hivatkozó sorok, pl. egy
% mentés, ugyanoda mutatnak), lefuttatja, és a végén (hiba esetén is) törli.
szoveg = fileread(bb84_fajl, 'Encoding', 'UTF-8');
for nev = fieldnames(par)'
    szoveg = regexprep(szoveg, param_minta(nev{1}), ['$1' mat2str(par.(nev{1}), 15) '$3'], ...
        'once', 'lineanchors');
end

[~, azon]  = fileparts(tempname);       % minden futásnál új név: a Matlab nem használhat régi változatot
ideiglenes = fullfile(fileparts(bb84_fajl), ['bb84_gui_' azon '.m']);
fid = fopen(ideiglenes, 'w', 'n', 'UTF-8');
if fid < 0
    error('bb84_gui:iras', 'Nem tudom létrehozni az ideiglenes fájlt: %s', ideiglenes);
end
fprintf(fid, '%s', szoveg);
fclose(fid);
torles = onCleanup(@() delete(ideiglenes));    % a függvény végén (hibánál is) törli
rehash;                                 % a Matlab lássa az új fájlt

[ered, naplo] = szkript_futtatasa(ideiglenes);
end


function [ered, naplo] = szkript_futtatasa(szkript_fajl) %#ok<INUSD> az evalc szövegében használja
% Külön függvény: a szkript (a clear-rel együtt) ennek a munkaterületén fut, így nem
% keveredik a felület változóival. A konzolra írt szöveget az evalc gyűjti össze.
naplo = evalc('run(szkript_fajl)');

% A szkript változói közül ezek kellenek a felületnek. Ami hiányzik, oda a szkript
% nem jutott el (pl. a QBER-küszöb miatt leállt), vagy a bb84.m nem számolja (pl.
% fajl_ered, ha nem ment).
kell = {'N_total', 'f_rep', 'tot_detected', 'tot_double', 'tot_sifted', 'tot_sifted_b', ...
    'QBER', 'QBER_b', 'a_m', 'fajl_ered', 'k_pe', 'Q_k', 'Q_up', 'Q_tol', ...
    'n_frames', 'n_kept', 'n_rec', 'leak_kept', 't_conf'};
ered = struct();
for iv = 1:numel(kell)
    if exist(kell{iv}, 'var') == 1
        ered.(kell{iv}) = eval(kell{iv});
    end
end
if isfield(ered, 'n_rec')
    % maradt-e hiba az egyeztetett kulcsban; a ka, kb a szkriptből jön, ezért kell az eval
    ered.hibas = eval('nnz(ka ~= kb)'); %#ok<EVLCS>
end
end


function eredmenyek(u, e, naplo, ido)
% Az eredmények megjelenítése: állapotsor, táblázat, kiírás
T_real = e.N_total / e.f_rep;
qs     = @(q, n) sprintf('%.3f ± %.3f', 100*q, 100*sqrt(q .* (1-q) ./ max(n, 1)));   % ± 1 szigma
szaz   = @(x) sprintf('%.3f', 100*x);
ap     = @(x) sprintf('%.3g', 100*x);      % az afterpulse-arány 1 MHz-en nagyon kicsi lehet

adat = {
    'Elküldött pulzusok',               ezres(e.N_total)
    'Detektált pulzusok',               ezres(e.tot_detected)
    'Dupla klikk (eldobva)',            ezres(e.tot_double)
    'Afterpulse / klikk [%]',           ap(e.a_m)
    'Siftelt bitek',                    ezres(e.tot_sifted)
    'Siftelt kulcsráta [bit/s]',        rata(e.tot_sifted / T_real)
    'QBER, Z-bázis [%]',                qs(e.QBER_b(1), e.tot_sifted_b(1))
    'QBER, X-bázis [%]',                qs(e.QBER_b(2), e.tot_sifted_b(2))
    'QBER, összesen [%]',               qs(e.QBER, e.tot_sifted)
    'Mintabitek (paraméterbecslés)',    ha_van(e, 'k_pe', @ezres)
    'Minta QBER, Q_k [%]',              ha_van(e, 'Q_k', szaz)
    'Serfling-korlát, Q_up [%]',        ha_van(e, 'Q_up', szaz)
    'Átment frame-ek',                  ha_van(e, 'n_kept', @(n) sprintf('%d / %d', n, e.n_frames))
    'Egyeztetett kulcs [bit]',          ha_van(e, 'n_rec', @ezres)
    'Egyeztetett kulcsráta [bit/s]',    ha_van(e, 'n_rec', @(n) rata(n / T_real))
    'Nyilvános bitek (paritás + hash)', ha_van(e, 'n_kept', @(n) ezres(e.leak_kept + e.t_conf * n))
    'Maradék hiba a kulcsban',          ha_van(e, 'hibas', @ezres)
    };
u.tabla.Data = adat;

% Állapotsor: meddig jutott el a szkript
if ~isfield(e, 'Q_up')
    uzenet = sprintf('Leállt: túl kevés a siftelt bit (%d). Növeld a szimulált időt.', e.tot_sifted);
elseif e.Q_k > e.Q_tol
    uzenet = sprintf('Leállt: a minta QBER-je (%.2f%%) a %.1f%%-os küszöb fölött van, nincs kulcs.', ...
        100*e.Q_k, 100*e.Q_tol);
elseif ~isfield(e, 'n_rec') || e.n_rec == 0
    uzenet = 'Leállt: egy frame sem ment át a megerősítésen, nincs kulcs.';
else
    uzenet = '';
end
if isempty(uzenet)
    u.allapot.Text = sprintf('Kész (%s): %s bites egyeztetett kulcs (%s bit/s)', ...
        ido_szoveg(ido), ezres(e.n_rec), rata(e.n_rec / T_real));
    if isfield(e, 'fajl_ered')          % ha a bb84.m elmenti az eredményeket
        [~, nev, kit] = fileparts(e.fajl_ered);
        u.allapot.Text = sprintf('%s, mentve: %s%s', u.allapot.Text, nev, kit);
    end
    u.allapot.FontColorMode = 'auto';
else
    u.allapot.Text = uzenet;
    u.allapot.FontColor = [208 59 59]/255;
end

u.naplo.Value = splitlines(strtrim(naplo));
end


function alap_visszaallit(fig)
% "Alapértékek a bb84.m-ből": a mezők visszaállítása (a bb84.m-et újra beolvassa)
u = fig.UserData;
try
    alap = alapertekek(u.bb84_fajl, u.P(:,1));
catch hiba
    uialert(fig, hiba.message, 'Nem sikerült beolvasni az alapértékeket');
    return;
end
for i = 1:size(u.P, 1)
    nev = u.P{i,1};
    ert = alap.(nev) / u.P{i,4};
    for d = 1:min(numel(ert), numel(u.mezo.(nev)))
        u.mezo.(nev)(d).Value = kerekit(ert(d));
    end
end
varhato_frissit(fig);
end


%% Apró segédfüggvények
function l = cimke(g, szoveg, sor, oszlop, igazitas)
l = uilabel(g, 'Text', szoveg, 'HorizontalAlignment', igazitas);
l.Layout.Row    = sor;
l.Layout.Column = oszlop;
end

function s = ha_van(e, nev, f)
% A szkript egy változója kiírva, vagy '–', ha a szkript oda nem jutott el
if isfield(e, nev)
    s = f(e.(nev));
else
    s = '–';
end
end

function x = kerekit(x)
% A mértékegység-váltás kerekítési hibája nélkül (pl. 50e-9 / 1e-9 = 50.000000000000007)
x = round(x, 10, 'significant');
end

function s = ezres(x)
% Egész szám ezres tagolással: 10000000 -> '10 000 000'
s = regexprep(sprintf('%d', round(x)), '(\d)(?=(\d{3})+$)', '$1 ');
end

function s = rata(x)
% Ráta kiírása: 13800 -> '13 800', 0.0412 -> '0.0412'
if x >= 1000
    s = ezres(x);
else
    s = sprintf('%.3g', x);
end
end

function s = ido_szoveg(t)
if t < 1
    s = '1 s alatt';
elseif t < 120
    s = sprintf('%.0f s', t);
elseif t < 7200
    s = sprintf('%.0f perc', t / 60);
else
    s = sprintf('%.1f óra', t / 3600);
end
end
