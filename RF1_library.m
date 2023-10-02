%RF1 library of funcitons

function [ZOCstubs]=Convert2OCStub(Z,Filtertype,StartElement)
    if Filtertype=="LP"
    %ZOCstubs reads as [Z OCStubZ0 OCStublength TLZ0 TLlength]
    %Works for ShuntCap with N=3
    if StartElement == "ShuntCap"
        [numRows,numCols] = size(Z);
        if numRows-1 == 3
            %Convert last OCStub and TL to TL and SCStub
            Identity="OCShuntStubTL2TLSCShuntStub";
            Z1=1;
            Z2=Z(numRows-1,3);
            [Z01,~,Z02,~]=KurodaTransform(Identity,Z1,Z2);
            %Convert last SCStub and TL to TL and OCStub
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z(numRows-2,2);
            Z2=Z01;
            [Z04,~,Z03,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z02;
            Z2=1;
            [Z06,~,Z05,~]=KurodaTransform(Identity,Z1,Z2);
            %Results
            ZOCstubs(1,1:2)=[Z(1,3) "OCShunt"];
            ZOCstubs(2,1:2)=[Z03 "TL"];
            ZOCstubs(3,1:2)=[Z04 "OCShunt"];
            ZOCstubs(4,1:2)=[Z05 "TL"];
            ZOCstubs(5,1:2)=[Z06 "OCShunt"];
        elseif numRows-1 == 4
            disp("N selection not supported.")
            ZOCstubs=0;
        elseif numRows-1 == 5
            %Step 1
            Identity="OCShuntStubTL2TLSCShuntStub";
            Z2=Z(5,3);
            Z1=1;
            [Z05,~,Z06,~]=KurodaTransform(Identity,Z1,Z2);
            %Step 2
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z(4,2);
            Z2=Z05;
            [Z05,~,Z04,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z06;
            Z2=1;
            [Z07,~,Z06,~]=KurodaTransform(Identity,Z1,Z2);
            %Step 3
            Identity="OCShuntStubTL2TLSCShuntStub";
            Z2=Z(3,3);
            Z1=Z04;
            [Z03,~,Z04,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="OCShuntStubTL2TLSCShuntStub";
            Z2=Z05;
            Z1=Z06;
            [Z05,~,Z06,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="OCShuntStubTL2TLSCShuntStub";
            Z2=Z07;
            Z1=1;
            [Z07,~,Z08,~]=KurodaTransform(Identity,Z1,Z2);
            %Step 4
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z(2,2);
            Z2=Z03;
            [Z03,~,Z02,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z04;
            Z2=Z05;
            [Z05,~,Z04,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z06;
            Z2=Z07;
            [Z07,~,Z06,~]=KurodaTransform(Identity,Z1,Z2);
            Identity="SCShuntStubTL2TLOCShuntStub";
            Z1=Z08;
            Z2=1;
            [Z09,~,Z08,~]=KurodaTransform(Identity,Z1,Z2);
            %Results
            ZOCstubs(1,1:2)=[Z(1,3) "OCShunt"];
            ZOCstubs(2,1:2)=[Z02 "TL"];
            ZOCstubs(3,1:2)=[Z03 "OCShunt"];
            ZOCstubs(4,1:2)=[Z04 "TL"];
            ZOCstubs(5,1:2)=[Z05 "OCShunt"];
            ZOCstubs(6,1:2)=[Z06 "TL"];
            ZOCstubs(7,1:2)=[Z07 "OCShunt"];
            ZOCstubs(8,1:2)=[Z08 "TL"];
            ZOCstubs(9,1:2)=[Z09 "OCShunt"];
        else
            disp("N selection not supported.")
            ZOCstubs=0;
        end
    elseif StartElement == "SeriesInd"

    else
        disp("error OCShuntConverter001");
    end
else
    disp("Filter selection not supported.")
    ZOCstubs=0;
end
end
function [Zshunts]=RLC2Zshunts(Z,Filtertype)
% this is a matrix conversion of the Richards transformation

if Filtertype == "LP"
    disp("Note: L becomes a SCstub and C becomes a OCstub");
    [numRows,numCols] = size(Z);
    x=1;
    while x <= numRows
        Zshunts(x,1)=Z(x,1);
        if Z(x,2)==0
            Zshunts(x,2)=Z(x,2);
        else
            Zshunts(x,2)=Z(x,2);
        end
        if Z(x,3)==0
            Zshunts(x,3)=Z(x,3);
        else
            Zshunts(x,3)=1/Z(x,3);
        end
        x=x+1;
    end
else
    disp("Filter selection not supported.")
    Zshunts=0;
end    
end
function [Z0,length]=RichardsTransform(LCval,Component)
syms lambda
% LCval=3;
% Component="Cap";    %"Cap" or "Ind"

if Component=="Cap"
    Z0=1/LCval;
    length=lambda/8;    %at Wc
elseif Component=="Ind"
    Z0=Ind;
    length=lambda/8;    %at Wc
else
    disp("error KurodaTransform001");
end

end
function [Z01,length1,Z02,length2]=KurodaTransform(Identity,Z1,Z2)
% Identity = "OCShuntStub2TL"; 
% "OCShuntStubTL2TLSCShuntStub" Open circuit shunt stub followed by a TL converted to TL followed by short circuit series stub
% "TLSCShuntStub2OCShuntStubTL" TL followed by short circuit series stub converted to Open circuit shunt stub followed by a TL
% "SCShuntStubTL2TLOCShuntStub" Short circuit series stub followed by a TL converted to TL followed by Open circuit shunt stub
% "TLOCShuntStub2SCShuntStubTL" TL followed by Open circuit shunt stub converted to Short circuit series stub followed by a TL

syms lambda
n2=1+(Z2/Z1);

if Identity == "OCShuntStubTL2TLSCShuntStub"
    length1=lambda/8;
    Z01=Z2/n2;
    length2=lambda/8;
    Z02=Z1/n2;
elseif Identity == "TLSCShuntStub2OCShuntStubTL"
    length1=lambda/8;
    Z01=Z2*n2;
    length2=lambda/8;
    Z02=Z1*n2;
elseif Identity == "SCShuntStubTL2TLOCShuntStub"
    length1=lambda/8;
    Z01=Z2*n2;
    length2=lambda/8;
    Z02=Z1*n2;
elseif Identity == "TLOCShuntStub2SCShuntStubTL"
    length1=lambda/8;
    Z01=Z1/n2;
    length2=lambda/8;
    Z02=Z1/n2;
else
    disp("error KurodaTransform001");
end
end
function [Znorm,Zdenorm]=TopicQuiz9LumpedHelp(f0,f2,f1,fa,ATTN,Filter,N,StartElement,Z0,Filtertype)
% fc=3e9;
% % BW=.1;
% % f1=f0-f0*BW/2
% % f2=f0+f0*BW/2
% f2=3.2e9;
% f1=2.8e9;
% fa=3.3e9;
% ATTN=30;
% Filter="EqualRipple3";   %"MaxFlat", "EqualRipple0.5", "EqualRipple3", and "MaxFlatTimeDelay"
% N=5;
% StartElement="ShuntCap";  %"ShuntCap" or "SeriesInd"
% Z0=50;
% Filtertype="BP";    %"LP", "HP", "BP", and "BS"
% [Znorm,Zdenorm]=TopicQuiz9LumpedHelp(fc,f2,f1,fa,ATTN,Filter,N,StartElement,Z0,Filtertype);
%
% Design an LC network for an RF bandpass filter with following characteristics:
%
% a) Center frequency: ω0 = 2π(3 GHz)
% b) ω2 = 2π(3.2 GHz), ω1 = 2π(2.8 GHz) and Bandwidth = ω2 − ω1 = 2π(0.4 GHz)
% c) Attenuation at ωa = 2π(3.3 GHz) is greater than 30 dB.
% d) Equal-ripple response with 3 dB ripples.
% e) Minimum number of elements (i.e. minimum order) for achieving the desired attenuation goal.
%
% In your solution, make sure to indicate the table and figure numbers used to design the low-pass prototype
% (LPP). Show a drawing of the LPP with normalized values for the elements. Use a network that starts
% with a shunt capacitor. Show a drawing of the transformed (in frequency, impedance, and bandpass
% response) LC network. For the LC network, make sure to calculate and indicate the inductor and capacitor
% values in pF and nH units.
%
% Hint: The order of the filter can be determined from low pass prototype curves by considering half of the
% bandpass response as if it is centered at DC. In other words, to determine the order of the filter, you can
% take ωc = ω2 − ω0 = 2π(0.2 GHz), ω = ωa − ω0 = 2π(0.3 GHz) and use the correct figure from
% figures 8.26 or 8.27 that corresponds to characteristics of your lowpass filter prototype (i.e. equal-ripple
% with 3 dB ripples).

W0=2*pi*f0;
W2=2*pi*f2;
W1=2*pi*f1;
Wa=2*pi*fa;
BW=W2-W1;
Wc=W2-W0;
W=Wa-W0;
W_Check=abs(W/Wc)-1;
disp("******************************************************************************");
String=sprintf("*** On the %s filter chart in the Pozar book pg 405 or 407",Filter);
disp(String);
String=sprintf("*** Confirm that line N=%d is at least %f dB at abs(W/Wc)-1=%f",N,ATTN,W_Check);
disp(String);
disp("******************************************************************************");

[G]=Gselect(N,Filter)
if Filtertype=="LP"    %"LP", "HP", "BP", and "BS"
    x=1;
    while x<=N+1
        [Orientation,ElementType,ElementTypeDetail,Component,Unit]=EvenOddElements(x,N,StartElement,Filtertype);
        ElementVal=G(1,x);
        if Component == "Cap"
            Znorm(x,1)=0;
            Znorm(x,2)=0;
            Znorm(x,3)=ElementVal;
            ElementVal=(ElementVal/(Z0*Wc))*1e12;
            Unit="pF";
            Zdenorm(x,1)=0;
            Zdenorm(x,2)=0;
            Zdenorm(x,3)=ElementVal;
        elseif Component == "Ind"
            Znorm(x,1)=0;
            Znorm(x,2)=ElementVal;
            Znorm(x,3)=0;
            ElementVal=(ElementVal*Z0/Wc)*1e9;
            Unit="nH";
            Zdenorm(x,1)=0;
            Zdenorm(x,2)=ElementVal;
            Zdenorm(x,3)=0;
        else
            Znorm(x,1)=ElementVal;
            Znorm(x,2)=0;
            Znorm(x,3)=0;
            ElementVal=ElementVal*Z0;
            Zdenorm(x,1)=ElementVal;
            Zdenorm(x,2)=0;
            Zdenorm(x,3)=0;
        end
        String=sprintf("Element %d is a %s %s with the value %f %s",x,Orientation,ElementType,ElementVal,Unit);
        disp(String);
        x=x+1;
    end
elseif Filtertype=="HP"
    Conversion=Filtertype;    %"LP", "HP", "BP", and "BS"
    x=1;
    while x<=N+1
        [Orientation,ElementType,ElementTypeDetail,Component,Unit]=EvenOddElements(x,N,StartElement,Filtertype);
        ElementVal=G(1,x);
        LCval=G(1,x);
        Znorm(x,1)=ElementVal;
        if Component == "Cap"
            [L,C]=ComponentConverter(Component,Conversion,LCval,Wc,W0,W1,W2);
            Znorm(x,1)=0;
            Znorm(x,2)=L;
            Znorm(x,3)=0;
            Ldenormalized=(L*Z0)*1e9;
            Zdenorm(x,1)=0;
            Zdenorm(x,2)=Ldenormalized;
            Zdenorm(x,3)=0;
            String=sprintf("Element %d is a %s %s %f nH",x,Orientation,ElementType,Ldenormalized);
            disp(String);
        elseif Component == "Ind"
            [L,C]=ComponentConverter(Component,Conversion,LCval,Wc,W0,W1,W2);
            Znorm(x,1)=0;
            Znorm(x,2)=0;
            Znorm(x,3)=C;
            Cdenormalized=(C/Z0)*1e12;
            Zdenorm(x,1)=0;
            Zdenorm(x,2)=0;
            Zdenorm(x,3)=Cdenormalized;
            String=sprintf("Element %d is a %s %s %f pF",x,Orientation,ElementType,Cdenormalized);
            disp(String);
        else
            Znorm(x,1)=ElementVal;
            Znorm(x,2)=0;
            Znorm(x,3)=0;
            Zdenormalized=ElementVal*Z0;
            Zdenorm(x,1)=Zdenormalized;
            Zdenorm(x,2)=0;
            Zdenorm(x,3)=0;
            String=sprintf("Element %d is a %s %s with the value %f %s",x,Orientation,ElementType,Zdenormalized,Unit);
            disp(String);
        end
        x=x+1;
    end
elseif Filtertype=="BP"
    Conversion=Filtertype;    %"LP", "HP", "BP", and "BS"
    x=1;
    while x<=N+1
        [Orientation,ElementType,ElementTypeDetail,Component,Unit]=EvenOddElements(x,N,StartElement,Filtertype);
        ElementVal=G(1,x);
        LCval=G(1,x);
        if Component == "Ind" || Component == "Cap"
            [L,C]=ComponentConverter(Component,Conversion,LCval,Wc,W0,W1,W2);
            Znorm(x,1)=0;
            Znorm(x,2)=L;
            Znorm(x,3)=C;
            Ldenormalized=(L*Z0)*1e9;
            Cdenormalized=(C/Z0)*1e12;
            Zdenorm(x,1)=0;
            Zdenorm(x,2)=Ldenormalized;
            Zdenorm(x,3)=Cdenormalized;
            String=sprintf("Element %d is a %s %s with a %f pF Cap and a %f nH Inductor",x,Orientation,ElementTypeDetail,Cdenormalized,Ldenormalized);
            disp(String);
        else
            Znorm(x,1)=ElementVal;
            Znorm(x,2)=0;
            Znorm(x,3)=0;
            Zdenormalized=ElementVal*Z0;
            Zdenorm(x,1)=Zdenormalized;
            Zdenorm(x,2)=0;
            Zdenorm(x,3)=0;
            String=sprintf("Element %d is a %s %s with the value %f %s",x,Orientation,ElementType,Zdenormalized,Unit);
            disp(String);
        end
        x=x+1;
    end
elseif Filtertype=="BS"
    disp("error FilterWork001");
else
    disp("error FilterWork002");
end

end
function [Orientation,ElementType,ElementTypeDetail,Component,Unit]=EvenOddElements(x,N,StartElement,Filtertype)
% x=1;
% N=5;
% StartElement="ShuntCap";  %"ShuntCap" or "SeriesInd"
% EvenOddElements(x,N,StartElement);
if Filtertype=="HP"
    if StartElement=="ShuntCap"
        if x==N+1
            Orientation="normalized";  %"Shunt" or "Series"
            ElementType="Impedance";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Skip";  %"Series LC", "Parallel LC", and "Skip"
            Component="Skip";    %"Ind", "Cap", and "Skip"
            Unit="Ohm"; %"F" or "H"
        elseif floor(x/2)==x/2
            % code for even
            Orientation="Series";  %"Shunt" or "Series"
            ElementType="Capacitor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Series LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Cap";    %"Ind", "Cap", and "Skip"
            Unit="H"; %"F" or "H"
        else
            % code for odd
            Orientation="Shunt";  %"Shunt" or "Series"
            ElementType="Inductor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Parallel LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Ind";    %"Ind", "Cap", and "Skip"
            Unit="F"; %"F" or "H"
        end
    elseif StartElement=="SeriesInd"
        if x==N+1
            Orientation="normalized";  %"Shunt" or "Series"
            ElementType="Impedance";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Skip";  %"Series LC", "Parallel LC", and "Skip"
            Component="Skip";    %"Ind", "Cap", and "Skip"
            Unit="Ohm"; %"F" or "H"
        elseif floor(x/2)==x/2
            % code for even
            Orientation="Shunt";  %"Shunt" or "Series"
            ElementType="Inductor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Parallel LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Ind";
            Unit="F"; %"F" or "H"
        else
            % code for odd
            Orientation="Series";  %"Shunt" or "Series"
            ElementType="Capacitor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Series LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Cap";    %"Ind", "Cap", and "Skip"
            Unit="H"; %"F" or "H"
        end
    end
else
    if StartElement=="ShuntCap"
        if x==N+1
            Orientation="normalized";  %"Shunt" or "Series"
            ElementType="Impedance";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Skip";  %"Series LC", "Parallel LC", and "Skip"
            Component="Skip";    %"Ind", "Cap", and "Skip"
            Unit="Ohm"; %"F" or "H"
        elseif floor(x/2)==x/2
            % code for even
            Orientation="Series";  %"Shunt" or "Series"
            ElementType="Inductor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Series LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Ind";    %"Ind", "Cap", and "Skip"
            Unit="H"; %"F" or "H"
        else
            % code for odd
            Orientation="Shunt";  %"Shunt" or "Series"
            ElementType="Capacitor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Parallel LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Cap";    %"Ind", "Cap", and "Skip"
            Unit="F"; %"F" or "H"
        end
    elseif StartElement=="SeriesInd"
        if x==N+1
            Orientation="normalized";  %"Shunt" or "Series"
            ElementType="Impedance";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Skip";  %"Series LC", "Parallel LC", and "Skip"
            Component="Skip";    %"Ind", "Cap", and "Skip"
            Unit="Ohm"; %"F" or "H"
        elseif floor(x/2)==x/2
            % code for even
            Orientation="Shunt";  %"Shunt" or "Series"
            ElementType="Capacitor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Parallel LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Cap";
            Unit="F"; %"F" or "H"
        else
            % code for odd
            Orientation="Series";  %"Shunt" or "Series"
            ElementType="Inductor";  %"Capacitor" or "Inductor"
            ElementTypeDetail="Series LC";  %"Series LC", "Parallel LC", and "Skip"
            Component="Ind";    %"Ind", "Cap", and "Skip"
            Unit="H"; %"F" or "H"
        end
    end
end

end
function [L,C]=ComponentConverter(Component,Conversion,LCval,Wc,W0,W1,W2)
% Component="Ind";    %"Ind" or "Cap"
% Conversion="BP";    %"LP", "HP", "BP", and "BS"
% LCval=3.4817;
% Wc=0.2e9;   %W2-W0
% W0=2*pi*3e9; %Center frequency
% W1=2*pi*2.8e9;  %Upper Cutoff
% W2=2*pi*3.2e9;  %Lower Cutoff
% [L,C]=ComponentConverter(Component,Conversion,LCval,Wc,W0,W1,W2)

delta=(W2-W1)/W0;
if Component=="Ind"
    if Conversion=="LP"
        L=LCval;
        C=0;
    elseif Conversion=="HP"
        L=0;
        C=1/(Wc*LCval);
    elseif Conversion=="BP"
        L=LCval/(W0*delta);
        C=delta/(W0*LCval);
    elseif Conversion=="BS"
        L=LCval*delta/W0;
        C=1/(W0*LCval*delta);
    else
        disp("error ComponentConverter001");
    end
elseif Component=="Cap"
    if Conversion=="LP"
        L=0;
        C=LCval;
    elseif Conversion=="HP"
        L=1/(Wc*LCval);
        C=0;
    elseif Conversion=="BP"
        L=delta/(W0*LCval);
        C=LCval/(W0*delta);
    elseif Conversion=="BS"
        L=1/(W0*LCval*delta);
        C=LCval*delta/W0;
        disp("error ComponentConverter002");
    end
else
    disp("error ComponentConverter003");
end
end
function [G]=Gselect(N,Filter)
%filters: "MaxFlat", "EqualRipple0.5", "EqualRipple3", and "MaxFlatTimeDelay"
if Filter=="MaxFlat"
    if N==1
        G = [2.0000 1.0000];
    elseif N==2
        G = [1.4142 1.4142 1.0000];
    elseif N==3
        G = [1.0000 2.0000 1.0000 1.0000];
    elseif N==4
        G = [0.7654 1.8478 1.8478 0.7654 1.0000];
    elseif N==5
        G = [0.6180 1.6180 2.0000 1.6180 0.6180 1.0000];
    elseif N==6
        G = [0.5176 1.4142 1.9318 1.9318 1.4142 0.5176 1.0000];
    elseif N==7
        G = [0.4450 1.2470 1.8019 2.0000 1.8019 1.2470 0.4450 1.0000];
    elseif N==8
        G = [0.3902 1.1111 1.6629 1.9615 1.9615 1.6629 1.1111 0.3902 1.0000];
    elseif N==9
        G = [0.3473 1.0000 1.5321 1.8794 2.0000 1.8794 1.5321 1.0000 0.3473 1.0000];
    elseif N==10
        G = [0.3129 0.9080 1.4142 1.7820 1.9754 1.9754 1.7820 1.4142 0.9080 0.3129 1.0000];
    else
        disp("error Gselect001");
    end
elseif Filter=="EqualRipple0.5"
    if N==1
        G = [0.6986 1.0000];
    elseif N==2
        G = [1.4029 0.7071 1.9841];
    elseif N==3
        G = [1.5963 1.0967 1.5963 1.0000];
    elseif N==4
        G = [1.6703 1.1926 2.3661 0.8419 1.9841];
    elseif N==5
        G = [1.7058 1.2296 2.5408 1.2296 1.7058 1.0000];
    elseif N==6
        G = [1.7254 1.2479 2.6064 1.3137 2.4758 0.8696 1.9841];
    elseif N==7
        G = [1.7372 1.2583 2.6381 1.3444 2.6381 1.2583 1.7372 1.0000];
    elseif N==8
        G = [1.7451 1.2647 2.6564 1.3590 2.6964 1.3389 2.5093 0.8796 1.9841];
    elseif N==9
        G = [1.7504 1.2690 2.6678 1.3673 2.7239 1.3673 2.6678 1.2690 1.7504 1.0000];
    elseif N==10
        G = [1.7543 1.2721 2.6754 1.3725 2.7392 1.3806 2.7231 1.3485 2.5239 0.8842 1.9841];
    else
        disp("error Gselect002");
    end
elseif Filter=="EqualRipple3"
    if N==1
        G = [1.9953 1.0000];
    elseif N==2
        G = [3.1013 0.5339 5.8095];
    elseif N==3
        G = [3.3487 0.7117 3.3487 1.0000];
    elseif N==4
        G = [3.4389 0.7483 4.3471 0.5920 5.8095];
    elseif N==5
        G = [3.4817 0.7618 4.5381 0.7618 3.4817 1.0000];
    elseif N==6
        G = [3.5045 0.7685 4.6061 0.7929 4.4641 0.6033 5.8095];
    elseif N==7
        G = [3.5182 0.7723 4.6386 0.8039 4.6386 0.7723 3.5182 1.0000];
    elseif N==8
        G = [3.5277 0.7745 4.6575 0.8089 4.6990 0.8018 4.4990 0.6073 5.8095];
    elseif N==9
        G = [3.5340 0.7760 4.6692 0.8118 4.7272 0.8118 4.6692 0.7760 3.5340 1.0000];
    elseif N==10
        G = [3.5384 0.7771 4.6768 0.8136 4.7425 0.8164 4.7260 0.8051 4.5142 0.6091 5.8095];
    else
        disp("error Gselect003");
    end
elseif Filter=="MaxFlatTimeDelay"
    if N==1
        G = [2.0000 1.0000];
    elseif N==2
        G = [1.5774 0.4226 1.0000];
    elseif N==3
        G = [1.2550 0.5528 0.1922 1.0000];
    elseif N==4
        G = [1.0598 0.5116 0.3181 0.1104 1.0000];
    elseif N==5
        G = [0.9303 0.4577 0.3312 0.2090 0.0718 1.0000];
    elseif N==6
        G = [0.8377 0.4116 0.3158 0.2364 0.1480 0.0505 1.0000];
    elseif N==7
        G = [0.7677 0.3744 0.2944 0.2378 0.1778 0.1104 0.0375 1.0000];
    elseif N==8
        G = [0.7125 0.3446 0.2735 0.2297 0.1867 0.1387 0.0855 0.0289 1.0000];
    elseif N==9
        G = [0.6678 0.3203 0.2547 0.2184 0.1859 0.1506 0.1111 0.0682 0.0230 1.0000];
    elseif N==10
        G = [0.6305 0.3002 0.2384 0.2066 0.1808 0.1539 0.1240 0.0911 0.0557 0.0187 1.0000];
    else
        disp("error Gselect004");
    end
else
        disp("error Gselect005");
    
end


end
function TopicQuiz8wilkinsonHelp(S,ZL,Z0)
syms Vpos
% S = [0 -1i/sqrt(2) -1i/sqrt(2);-1i/sqrt(2) 0 0;-1i/sqrt(2) 0 0]
% ZL=100;
% Z0=50;
% TopicQuiz8wilkinsonHelp(S,ZL,Z0)
%
% Topic Quiz 8 Question
% S-parameters of Wilkinson power divider is referenced to 50Ω port impedance and given as:
%
% S = [0 -1i/sqrt(2) -1i/sqrt(2);-1i/sqrt(2) 0 0;-1i/sqrt(2) 0 0]
%
% Port 1 is terminated with a 100 Ω load. Incident voltages at port 2 and port 3 are given as V2+ = V3+ = 1V.
% Determine power delivered to port 1 and reflected powers from port 2 and port 3.

[numRows,numCols] = size(S);
R=1;
C=1;
Vneg=zeros(numRows,1);

while R <= numRows
    while C <= numCols
        Vneg(R,1)=Vneg(R,1)+S(R,C);
        C=C+1;
    end
    C=1;
    R=R+1;
end

RefCoeff=(ZL-Z0)/(ZL+Z0);
R=1;
while R <= numRows
    if R==1
%         Vpos=solve(RefCoeff==Vpos/Vneg(R,1),Vpos);
        Vpos=Vneg(R,1)*RefCoeff;
        Vtot=Vpos+Vneg(R,1);
        Vtot_rounded=round(Vtot,4);
        String=sprintf("Vtot at port %d is",R);
%         disp(String);
%         disp(Vtot_rounded);
        Itot=-(Vpos-Vneg(R,1))/Z0;
        Itot_rounded=round(Itot,4);
        String=sprintf("Itot at port %d is",R);
%         disp(String);
%         disp(Itot_rounded);
        PL=0.5*real(Vtot*conj(Itot));
        String=sprintf("Power reflected at port %d is %f",R,PL);
        disp(String);
        disp(" ");
    else
%         disp("row "+R);
        Vneg_otherrows=Vneg(R,1)*Vpos;
        String=sprintf("V- at port %d is %f",R,Vneg_otherrows);
%         disp(String);
        PL=0.5*real(Vneg_otherrows*conj(Vneg_otherrows/Z0));
        String=sprintf("Power reflected at port %d is %f",R,PL);
        disp(String);
        disp(" ");
    end
    R=R+1;
end

end
function TopicQuiz7matchHelp(ZL,Zin,start)
% ZL=25-50i;  
% Zin=50;
% start="SeriesInd"; %"SeriesInd","SeriesCap","ShuntInd", and "ShuntCap"
% TopicQuiz7matchHelp(ZL,Zin,start);

%Topic Quiz 7 Question
% Part (a) Design an L-section matching network to transform ZL = 25 − j50 Ω to Zin = 50 Ω.
%
% Your matching network should start with a series inductor. If there are multiple choices for a series
% inductor, use the inductor that will require smallest reactance. Draw the matching network’s circuit
% schematic starting with the load on the right hand side of the network. Your drawing should indicate the
% type of series and shunt elements (capacitors or inductors). Calculate the un-normalized reactance for
% both series and shunt element in your matching network.
%
% Solution should be shown on the type of Smith Chart provided below.
%
% Part (b) For the same matching network problem, how many different L-section matching network
% designs are possible? For this part, use a new Smith Chart and identify all matching network possibilities.
% You do not need to calculate the reactance for the elements in your matching network. Only show the
% circuit diagram for each matching network – your drawing should show/indicate if circuit elements in the
% network are capacitors or inductors.

Znorm=Normalize(ZL,Zin);
if start=="SeriesInd"
    disp("up right");
elseif start=="SeriesCap"
    disp("down right");
elseif start=="ShuntInd"
    disp("up left");
elseif start=="ShuntCap"
    disp("down left");
else
    disp("error");
end

disp("remember if using the shunt movement on the smith chart you need to put your answer under 1 to get impedance");

end
function [RefCoeffmag,RefCoeffAng]=TopicQuiz6SmithHelp(ZTL,length,Z,dir,mag,angle)
syms lambda
%Topic Quiz 6 Question
%Find the Γin and Zin using the Smith Chart.
% ZTL=100;          %characteristic impedance of TL
% length=lambda/3;  %lambda must be syms
% Z=50+50i;         %impedance of load or generator
% dir="gen";        %"gen" or "load"
% mag=1.5;          %inches
% angle=116;        %need to calc from smithshart
% [RefCoeffmag,RefCoeffAng]=TopicQuiz6SmithHelp(ZTL,length,Z,dir,mag,angle)

disp("Step01: normalize impedances");
Znorm=Normalize(Z,ZTL);
disp("Step02: Plot on smith chart");
disp("Step03: Use compass to draw circle with constant reflection coefficient (magnitude)");
disp("Step04: draw line from center through impedance value");
disp("Step05: rotate clockwise if moving toward generator and counter clockwise otherwise");
if dir=="gen"
    Value=-1;
    disp("move clockwise");
else
    Value=1;
    disp("move counter clockwise");
end
disp("Step06: remember beta*l=(2pi/lambda)*l");
Angrad=(2*pi/lambda)*length;
Angdeg=Angrad*(180/pi);
disp("Step07: double the above value");
Ang=Angdeg*2;
disp("Step08: rotate the above value in the direction determined in step 5");
RefCoeffAng=round(angle+Value*Ang,4)
disp("Step09: This is the new value in the form mag(Z)e^jangle(Z)");
RefCoeffmag=mag*1/3.25
end
function Znorm=Normalize(Z,Z0)
%Z=25+i50;
%Z0=50;
%Znorm=Normalize(Z,Z0);
%Basic function
%Z0 is usually 50 Ohms

Znorm=Z/Z0;


end
function N5_3db_Filter(fc,percentage,ATTN_Freq,Z0) %Cascading
% fc=3e9
% percentage=0.035 %0.023101     %0 to 2.3101%
% ATTN_Freq=3.55e9
% Z0=50;
% N5_3db_Filter(fc,percentage,ATTN_Freq,Z0);

w0=2*pi*fc;
w1=w0*(1-percentage/2);
f1=w0*(1-percentage/2)/(2*pi)
w2=w0*(1+percentage/2);
f2=w0*(1+percentage/2)/(2*pi)
del=(w2-w1)/w0;
wa=2*pi*ATTN_Freq;
w=wa-w0;
wc=w2-w0;
Norm_freq=abs(w/wc)-1;

%N=5 Element values for 3 dB ripple
g1=3.4817;
g2=0.7618;
g3=4.5381;
g4=0.7618;
g5=3.4817;
g6=1.0000;

J1Z0=sqrt(pi*del/(2*g1));
J2Z0=pi*del/(2*sqrt(g1*g2));
J3Z0=pi*del/(2*sqrt(g2*g3));

Z0e1=Z0*(1+J1Z0+J1Z0^2)
Z0e2=Z0*(1+J2Z0+J2Z0^2)
Z0e3=Z0*(1+J3Z0+J3Z0^2)
Z0o1=Z0*(1-J1Z0+J1Z0^2)
Z0o2=Z0*(1-J2Z0+J2Z0^2)
Z0o3=Z0*(1-J3Z0+J3Z0^2)
end
function N4_MAXFLAT_LPFilter(fc,ATTN_Freq,Z0)
% fc=2e9;
% ATTN_Freq=3.4e9;
% Z0=50;
% N4_MAXFLAT_LPFilter(fc,ATTN_Freq,Z0);

Wc=2*pi*fc;
Wa=2*pi*ATTN_Freq;
Norm_freq=abs(Wa/Wc)-1

%N=5 Element values for Max Flat
g1=0.7654;
z1=g1;
g2=1.8478;
z2=1/g2;
g3=1.8478;
z3=g3;
g4=0.7654;
z4=1/g4;
g5=1;

N=1+1/z1;
Z1=N*z1
Z2=1/(N*z2)
N=1+z4/1;
Z3=z3/N
Z4=z4/N
N1=1+Z4/g3
N2=1+1/Z3

% C1_prime=g1/(Z0*Wc)
% L2_prime=Z0*g2/Wc
% C3_prime=g3/(Z0*Wc)
% L4_prime=Z0*g4/Wc
end
function N5_Equal_HPFilter(fc,ATTN_Freq,Z0)
% fc=3e9;
% ATTN_Freq=2e9;
% Z0=75;
% N5_Equal_HPFilter(fc,ATTN_Freq,Z0);

Wc=2*pi*fc;
Wa=2*pi*ATTN_Freq;
Norm_freq=abs(Wc/Wa)-1

%N=5 Element values for Max Flat
g1=3.4817;
g2=0.7618;
g3=4.5381;
g4=g2;
g5=g1;
g6=1.0000;

C1_prime=Z0/(g1*Wc)
L2_prime=1/(Z0*g2*Wc)
C3_prime=Z0/(g3*Wc)
L4_prime=1/(Z0*g4*Wc)
C5_prime=Z0/(g5*Wc)

end
function N5_MAXFLAT_LPFilter(fc,ATTN_Freq,Z0)
% fc=2e9;
% ATTN_Freq=3.4e9;
% Z0=50;
% N5_MAXFLAT_LPFilter(fc,ATTN_Freq,Z0);

Wc=2*pi*fc;
Wa=2*pi*ATTN_Freq;
Norm_freq=abs(Wa/Wc)-1

%N=5 Element values for Max Flat
g1=0.6180;
g2=1.6180;
g3=2.0000;
g4=1.6180;
g5=0.6180;
g6=1.0000;

C1_prime=g1/(Z0*Wc)
L2_prime=Z0*g2/Wc
C3_prime=g3/(Z0*Wc)
L4_prime=Z0*g4/Wc
C5_prime=g5/(Z0*Wc)

end
function DirCoupSolver(P1_dBm,C,D,L)
% P1_dBm=20 %dBm
% C=20;   %dB
% D=35;   %dB
% L=0.5;  %dB

P1_mW=10^(P1_dBm/10)   %mW
P2_mW=P1_mW/(10^(L/10))
P2_dBm=10*log10(P2_mW/1)
P3_mW=P1_mW/(10^(C/10))
P3_dBm=10*log10(P3_mW/1)
P4_mW=P3_mW/(10^(D/10))
P4_dBm=10*log10(P4_mW/1)
end
function N3_3db_Filter(fc,percentage,ATTN_Freq,Z0) %Cascading
% fc=3e9
% percentage=0.023101     %0 to 2.3101%
% ATTN_Freq=3.55e9
% Z0=50;

w0=2*pi*fc;
w1=w0*(1-percentage/2);
f1=w0*(1-percentage/2)/(2*pi);
w2=w0*(1+percentage/2);
f2=w0*(1+percentage/2)/(2*pi);
del=(w2-w1)/w0
wa=2*pi*ATTN_Freq;
w=wa-w0;
wc=w2-w0;
Norm_freq=abs(w/wc)-1

%N=3 Element values for 3 dB ripple
g1=3.3487;
g2=0.7117;
g3=3.3487;
g4=1;

% g1=1.5963;
% g2=1.0967;
% g3=1.5963;
% g4=1;

L1=del*Z0/(w0*g1);
C1=g1/(del*w0*Z0);
L2=g2*Z0/(del*w0);
C2=del/(w0*g2*Z0);

J1Z0=sqrt(pi*del/(2*g1));
J2Z0=pi*del/(2*sqrt(g1*g2));
% J3Z0=pi*del/(2*sqrt(g2*g3));

Z0e1=Z0*(1+J1Z0+J1Z0^2)
Z0e2=Z0*(1+J2Z0+J2Z0^2)
Z0o1=Z0*(1-J1Z0+J1Z0^2)
Z0o2=Z0*(1-J2Z0+J2Z0^2)
end
function [Q, R] = Shmitty(A) %Cascading
% A=[1 0 0; 0 1 0; 0 0 1]
a=A(:,1);
b=A(:,2);
c=A(:,3);
 
A = a;
B = b-(A'*b/(A'*A))*A;
C = c-(A'*c/(A'*A))*A-(B'*c/(B'*B))*B;
q1 = A./norm(A);
q2 = B./norm(B);
q3 = C./norm(C);
Q = [q1 q2 q3]
R = [q1'*a q1'*b q1'*c; 0 q2'*b q2'*c; 0 0 q3'*c]
end
function [A,B,C,D,ABCD] = S2ABCD(S11,S12,S21,S22,Z0) %Cascading
% S11 = 1;
% S12 = 2;
% S21 = 3;
% S22 = 4;
% Z0 = 50;
% [A,B,C,D,ABCD] = S2ABCD(S11,S12,S21,S22,Z0)
A = ((1+S11)*(1-S22)+S12*S21)/(2*S21);
B = Z0*((1+S11)*(1+S22)-S12*S21)/(2*S21);
C = (1/Z0)*((1-S11)*(1-S22)-S12*S21)/(2*S21);
D = ((1-S11)*(1+S22)+S12*S21)/(2*S21);
ABCD = [A B;
        C D];
% disp(ABCD);
end
function [S11,S12,S21,S22,S] = ABCD2S(A,B,C,D,Z0) %Cascading
% A = 1;
% B = 2;
% C = 3;
% D = 4;
% Z0 = 50;
% [S11,S12,S21,S22,S] = ABCD2S(A,B,C,D,Z0)
S11 = (A+B/Z0-C*Z0-D)/(A+B/Z0+C*Z0+D);
S12 = (2*(A*D-B*C))/(A+B/Z0+C*Z0+D);
S21 = (2)/(A+B/Z0+C*Z0+D);
S22 = (-A+B/Z0-C*Z0+D)/(A+B/Z0+C*Z0+D);
S = [S11 S12;
     S21 S22];
% disp(S);
end
function [A,Theta] = Rect2Pol(R,I)
% R = 35;   %Real R
% I = 16;  %imaginary X
A = sqrt((R.^2)+(I.^2));
Theta = atan(I./R);
if R<0 && Theta>=-pi
    Theta=Theta+pi;
elseif R<0 && Theta<=pi
    Theta=Theta-pi;
end
% Thetad = atand(I/R);
% disp("Polar value is "+A+" at angle "+Theta+" rad or "+Thetad+" deg");
end
function [A,Theta] = Rect2Pold(R,I)
% R = 35;   %Real R
% I = 16;  %imaginary X
A = sqrt((R.^2)+(I.^2));
Theta = atand(I./R);
if R<0 && Theta>=-180
    Theta=Theta+180;
elseif R<0 && Theta<=180
    Theta=Theta-180;
end
% Thetad = atand(I/R);
% disp("Polar value is "+A+" at angle "+Theta+" rad or "+Thetad+" deg");
end
function Rect = Pol2Rect(A,Theta)
% A = 35;   %Amplitude
% Theta = pi/3;  %rad
R = A*cos(Theta);
I = A*sin(Theta);
Rect = R + I*1i;
% disp("Rectangular value is "+Rect);
end
function Rect = Pol2Rectd(A,Theta)
% A = 35;   %Amplitude
% Theta = pi/3;  %rad
R = A*cosd(Theta);
I = A*sind(Theta);
Rect = R + I*1i;
% disp("Rectangular value is "+Rect);
end
function CharImpedance = CharImp(L,R,w,PropConstant)
% L = 35e-6;   %H/m
% R = 4;  %Ohm/m
% w = 2*pi*(800e6); %rad/sec
% PropConstant = 10+10i;    %units
CharImpedance = (R+1i*w*L)/PropConstant;
CharImpedanceR = round(CharImpedance,2,'decimals');
disp("Characteristic impedance is "+CharImpedanceR);
end
function PropConstant = PropCons(L,C,R,G,w)
% L = 35e-6;   %H/m
% C = 25e-12;   %F/m
% R = 4;  %Ohm/m
% G = 0.02; %S/m
% w = 2*pi*(800e6); %rad/sec
PropConstant = sqrt((R+1i*w*L)*(G+1i*w*C));
PropConstantR = round(PropConstant,2,'decimals');
disp("Propigation constant is "+PropConstantR);
end
function Zin1 = ZinCalcAll(Z0,ZL,Beta,L)   %Not Quarter Wavelength
% Z0 = 35;
% ZL = 25 + 1i*25;
% Beta = 2*pi/lambda;
% L = lambda/8;
Zin1 = Z0*((ZL + 1i*Z0*tan(Beta*L))/(Z0 + 1i*ZL*tan(Beta*L)));
Zin1R = round(Zin1,2,'decimals');
disp(Zin1R);
end
function Zin2 = ZinCalcQrt(Z0,ZL)   %Quarter Wavelength
% Z0 = 65;
% ZL = Zin1;
% Beta = 2*pi/lambda;
% L = lambda/4;
Zin2 = (Z0^2)/ZL;
Zin2R = round(Zin2,2,'decimals');
disp(Zin2R);
end
function Zout = VDiv(V,R1,R2)
% V = 10;
% R1 = 80-1i*40;
% R2 = 100;
Zout = V*(R2)/(R1+R2);
ZoutR = round(Zout,2,'decimals');
disp(ZoutR);
end