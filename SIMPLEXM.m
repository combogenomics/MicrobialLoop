function [x,options] = SIMPLEXM(funfcn,x,options,grad,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10);
%SIMPLEX	versione 1.5
%Trova il minimo di una funzione a pi? variabili utilizzando il metodo del simplesso
% con la fase di espansione ottimizzata.
%
%Sintassi minima: X=SIMPLEX('FUN',X0) parte dalle condizioni iniziali X0 e 
% trova il minimo della funzione descritta in FUN (usualmente un M-file di nome FUN.M).
%	La funzione FUN deve ritornare un valore scalare.
%
%Sintassi completa: [X,OPTIONS]=SIMPLEX('FUN',X0,OPTIONS,[],P1,P2,..,P10)
%	OPTIONS: ? un vettore di parametri opzionali (vedere FOPTIONS()).
%	OPTIONS(1): permette di visualizzare i risultati intermedi
%		0 non visualizza i risultati intermedi
%		1 visualizza i risultati intermedi
%		2 visualizza i risultati intermedi in forma ridotta
%		  (numero di iterazioni e valore del funzionale di errore)
%	OPTIONS(2) ? la precisione richiesta per il valore di X alla soluzione.
%	OPTIONS(3) ? la precisione richiesta per il valore del funzionale di errore
%		alla soluzione.
%	OPTIONS(6): considera il file STARTV per i vertici iniziali del simplesso
%		0 costruisce il simplesso di partenza da X0.
%		1 reinizializza i vertici del simplesso dai dati memorizzati in STARTV
%		2 terminata l'ottimizzazione salva i vertici del simplesso in STARTV
%		3 utilizza STARTV sia per l'inizializzazione del simplesso iniziale
%		  sia per salvare il simplesso finale.
%	OPTIONS(14) ? il numero massimo di iterazioni dopo il quale la ricerca
%		del minimo si blocca. Se ? negativa indica dopo quanti minuti la ricerca
%		si deve bloccare.
%		alla soluzione.
%	P1,P2,....,P10 sono delle variabili da passare direttamente a FUN().
x
funfcn
%P1=P1;P2=P2;P3=P3;
if nargin<3, options=[]; end
options=foptions(options);
prnt=options(1);
tol=options(2);
tol2=options(3);
evalstr = [funfcn];
if ~any(funfcn<48)
	evalstr=[evalstr, '(x'];
	for i=1:nargin - 4
		evalstr = [evalstr,',P',num2str(i)];
	end
	evalstr = [evalstr, ')'];
end

n = prod(size(x));
if (~options(14))
	options(14) = 200*n; 
end

% Inizializza il simplesso
if (options(6)==1 || options(6)==3)
	if prnt
		disp('Ripristino del simplesso di partenza')
	end
	%Riprendi il simplesso memorizzato
	load STARTV;
else
	%Inizializza il simplesso partendo dalla condizione iniziale 
	if prnt
		disp('Calcolo del simplesso di partenza')
	end
	xin = x(:)
	v = 0.8*xin
	x(:) = v
	 fv = eval(evalstr)
	for j = 1:n
		y = xin;
		if y(j) ~= 0
			y(j) = 1.2*y(j);
		else
			y(j) = 0.1;
		end
		v = [v y];
		x(:) = y; f = eval(evalstr);
		fv = [fv  f];
	end
	[fv,j] = sort(fv);
	v = v(:,j);
end
cnt = n+1;

% nr=1000;
% for iir=1:nr,
%     vnew(:,iir)=v(:,1)+0.25*v(:,1)*randn(1);
%     x(:)=vnew(:,iir);
%     fnew=eval(evalstr);
%     FNEW(:,iir)=fnew;
% end
% [fvnr,jjr]=sort(FNEW);
% vnr=vnew(:,jjr);
% v=vnr(:,1:n+1);fv=fvnr(1:n+1);
% cnt=cnt+nr;

if prnt
	clc
	format compact
	format short e
	home
	cnt
	disp('Iniziale')
	disp(' ')
	v
	fv
end

alpha=1.5;beta=0.5;gamma=2;fib_max=0;
%alpha=0.8;beta=0.4;gamma=1.4;fib_max=0;
%alpha=1.25;beta=0.5;gamma=1.75;fib_max=0;
contatori=[0 0 0 0 0];
[n,np1] = size(v);
onesn = ones(1,n); 
ot = 2:n+1;
on = 1:n;
tStart=clock;
fOverTime=0;
fOverCnt=0;
fContinua=1;

% alfa=20;conta=0;
% VRINT=zeros(n,alfa*(n+1));
% FRINT=zeros(1,alfa*(n+1));
% while alfa>=1,
% conta=conta+1;
% alfa=alfa-1;   
% vbar=(sum(v(:,on)')/n)'; % calcolo baricentro
% [vv,fvv]=multidirectional(x,v,fv,vbar,evalstr,n,on,alfa,P1,P2,P3);cnt=cnt+1;
% if conta==1,
% VRINT(:,1:n+1)=vv;
% FRINT(1:n+1)=fvv;
% else
% VRINT(:,(conta-1)*(n+1)+1:conta*(n+1))=vv;
% FRINT((conta-1)*(n+1)+1:conta*(n+1))=fvv;
% end
% end
% [SFRINT sk]=sort([FRINT]);
% Vstot=VRINT(:,sk);
% Vmstot=Vstot(:,1:n+1);
% SFNEW=SFRINT(1:n+1);
% v=Vmstot;
% fv=SFNEW;
PS=10;
%Loop principale.
%while fContinua & v>-100 & v<100,%cnt < options(14)
while fContinua %& v>-50 & v<50, 
	vbar=(sum(v(:,on)')/n)'; % calcolo nuovo baricentro
%     %Riflessione primaria
% 	vr=vbar+alpha*(vbar-v(:,n+1)); x(:) = vr; fr = eval(evalstr);cnt=cnt+1;%riflessione tramite il baricentro
% 	vk = vr; fk = fr; how = 'RIFLESSIONE';azione=1;
	%Riflessione primaria multipla (si valutano più riflessioni e si prende la migliore)
    VRINT=zeros(n,PS);
    FRINT=zeros(1,PS);
    for iir=1:PS,
	     %vr=vbar+((alpha-0.4)+0.4*(iir/PS))*(vbar-v(:,n+1));
         vr=vbar+((alpha-1.4)+1.4*(iir/PS))*(vbar-v(:,n+1));
         x(:) = vr; fr = eval(evalstr);cnt=cnt+1;%riflessione tramite il baricentro
         VRINT(:,iir)=vr;
         FRINT(iir)=fr;
    end
    [SFRINT sk]=sort(FRINT);
    Vtot=VRINT;
    Vstot=Vtot(:,sk);
    Vmstot=Vstot(:,1);
    SFNEW=SFRINT(1);
    vr=Vmstot;
    fr=SFNEW;
	vk = vr; fk = fr; how = 'RIFLESSIONE';azione=1;
	if fr<fv(n) && fr~=NaN
		if fr<fv(1)
            %vr2=vbar+gamma*(vr-vbar);x(:)=vr2;fr2=eval(evalstr);cnt=cnt+1;how='ESPANSIONE';azione=2;
			%Espansione multipla (si valutano più espansioni e si prende la migliore)
            VSINT=zeros(n,PS);
            FSINT=zeros(1,PS);
            for iie=1:PS,
			     %vr2=vbar+((gamma-0.5)+0.5*(iie/PS))*(vr-vbar);
                 vr2=vbar+((gamma-1)+1*(iie/PS))*(vr-vbar);
                 x(:)=vr2;fr2=eval(evalstr);cnt=cnt+1;how='ESPANSIONE';azione=2;
                 VSINT(:,iie)=vr2;
                 FSINT(iie)=fr2;
            end
            [SFSINT sks]=sort(FSINT);
            Vetot=VSINT;
            Vestot=Vetot(:,sks);
            Vsstot=Vestot(:,1);
            SSNEW=SFSINT(1);
            vr2=Vsstot;
            fr2=SSNEW;
			if fr2<fv(1) && fr2~=NaN %Espansione ottimizzata
				search=1;rapp=1.5;ifib=2;
				while search %cerca intervallo
					vr3=vr2+rapp*(vr2-vr);x(:)=vr3;fr3=eval(evalstr);cnt=cnt+1;
					if fr3>fr2 && fr3~=NaN
						search=0;
					else
						vr=vr2;fr=fr2;
						vr2=vr3;fr2=fr3;
						ifib=ifib+1;
					end
				end
				%Il minimo ? compreso fra vr e vr3 => riduzione dell'intervallo
				if ifib>2
					%Calcolo vr2b in modo che lungo la direzione di ricerca si trovino nell'ordine vr,vr2,vr2b,vr3.
					vr2b=vr+rapp*(vr2-vr);x(:)=vr2b;fr2b=eval(evalstr);cnt=cnt+1;
					rapp=0.5;
					for fr=ifib:-1:3
						if fr2b<fr2
							vr=vr2;
							vr2=vr2b;fr2=fr2b;
							vr2b=vr3-rapp*(vr3-vr);x(:)=vr2b;fr2b=eval(evalstr);
						else
							vr3=vr2b;
							vr2b=vr2;fr2b=fr2;
							vr2=vr+rapp*(vr3-vr);x(:)=vr2;fr2=eval(evalstr);
						end
						cnt=cnt+1;
					end
					if fr2b<fr2
						vr2=vr2b;fr2=fr2b;
					end
				end
				vk=vr2;fk=fr2;how='ESPANSIONE OTTIMIZZATA';azione=3;
			end
		end
	else
		if max(max(abs(v(:,ot)-v(:,onesn)))) <= tol &&  max(abs(fv(1)-fv(ot))), break, end
%         if max(max(abs(v(:,ot)-v(:,onesn)))) <= tol &  max(abs(fv(1)-fv(ot))) <= tol2...
%                 & (x(1)>5 & x(1)<7 & x(2)>20 & x(2)<25 & x(3)>0.5 & x(3)<1 & x(4)>0 & x(4)<0.001 & x(5)>0.5 & x(5)<1 &...
%                 x(6)>2.3 & x(6)<2.8 & x(7)>0.25 & x(7)<0.6), break, end
		vt = v(:,n+1); ft = fv(n+1);  % vedi riga 115
		if fr < ft
			vt = vr;  ft = fr;        % vedi riga 115
		end
		vc=vbar+beta*(vt-vbar);x(:)= vc;fc=eval(evalstr);cnt=cnt+1;
		if fc<fv(n+1)
			vk = vc; fk = fc;
			how = 'CONTRAZIONE';azione=4;
		else
			for j = 2:n
				v(:,j)=(v(:,1)+v(:,j))/2;x(:)=v(:,j);fv(j)=eval(evalstr);
			end
			cnt = cnt + n-1;
			vk = (v(:,1) + v(:,n+1))/2; x(:) = vk; fk = eval(evalstr);cnt=cnt+1;
			how = 'RIDUZIONE';azione=5;
		end;
	end

	v(:,n+1)=vk;
	fv(n+1)=fk;
	[fv,j]=sort(fv);
	v=v(:,j);
	contatori(azione)=contatori(azione)+1;

	if prnt==1
		home;clc
		cnt
		disp(how)
		v
		fv
	elseif prnt==2
		home;clc
		fprintf('Contatore: %g           Funzionale di errore minimo: %g\n',cnt,fv(1));
	end;
	%Controlla la condizione di file loop
	if options(14)>=0
		if cnt > options(14)
			fContinua=0;
			fOverCnt=1;
		end
	else
		if etime(clock,tStart)>-options(14)*60
			fContinua=0;
			fOverTime=1;
		end
	end
end
x(:) = v(:,1);
if prnt, format, end
options(10)=cnt;
options(8)=min(fv); 
if (options(6)==2 || options(6)==3)
	%Salva il simplesso
	save STARTV v fv;
end

if fOverCnt
	if options(1) >= 0
		disp(['Attenzione: ? stato raggiunto il massimo numero di iterazioni previste (',num2str(options(14)),')'])
		disp( '         (Aumentare OPTIONS(14)).')
	end
elseif fOverTime
	if options(1) >= 0
		disp(['Attenzione: ? stato superato il limite massimo di tempo (',num2str(-options(14)),' minuti)'])
		disp( '         (Diminuire OPTIONS(14)).')
	end
end

if prnt
	fprintf('Numero riflessioni : %g\n',contatori(1))
	fprintf('Numero espansioni : %g\n',contatori(2))
	fprintf('Numero espansioni ottimizzate : %g\n',contatori(3))
	fprintf('Numero contrazioni : %g\n',contatori(4))
	fprintf('Numero riduzioni : %g\n',contatori(5))
end

