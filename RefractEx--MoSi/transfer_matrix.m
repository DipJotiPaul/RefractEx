% Initialization: thickness=[NbN Silicon], polarization: 0 for TM (-p) and 90 for TE (-s)
% for air-NbN-Silicon-air, layers=4, nk_layer=values(4,len of wv)
function [Tran,Ref,Abs]=transfer_matrix(wv,incidence_angle,polarization,thickness,nk_layer)
    layers=length(thickness)+2;
    if(max(thickness)<1)                                 % coherent layers formalism when thickness is below 1um
        if(polarization==0 || polarization==90) 
            [Tran,Ref]=tmm(wv,incidence_angle,polarization,thickness,nk_layer,layers);
        else
            [Tpol,Rpol]=tmm(wv,incidence_angle,0,thickness,nk_layer,layers);
            [Tsol,Rsol]=tmm(wv,incidence_angle,90,thickness,nk_layer,layers);
            Tran=0.5*(Tpol+Tsol);
            Ref=0.5*(Rpol+Rsol);
        end
    else                                                               % incoherent layers formalism in random phase when thickness > 1um
        Tran=zeros(1,length(wv));       Ref=zeros(1,length(wv));
        iteration=500;
        for k1=1:iteration
            if(polarization==0 || polarization==90)
                [T,R]=tmm(wv,incidence_angle,polarization,thickness,nk_layer,layers);
            else
                [Tpol,Rpol]=tmm(wv,incidence_angle,0,thickness,nk_layer,layers);
                [Tsol,Rsol]=tmm(wv,incidence_angle,90,thickness,nk_layer,layers);
                T=0.5*(Tpol+Tsol);
                R=0.5*(Rpol+Rsol);
            end
        Tran=Tran+T;        Ref=Ref+R;    
        end
        Tran=smoothdata(Tran/iteration,'gaussian',70);   
        Ref=smoothdata(Ref/iteration,'gaussian',70);    
    end
    Tran=Tran*100;      Ref=Ref*100;        Abs=100-Tran-Ref;
end

function [Tran,Ref]=tmm(wv,incidence_angle,polarization,thickness,nk_layer,layers)
    Tran=zeros(1,length(wv));       Ref=zeros(1,length(wv));
    for k1=1:length(wv)
        angle=[incidence_angle zeros(1,layers-1)];
        for k2=1:layers-1
            [coherent_matrix,angle(k2+1)]=coherent_tmm(wv(k1),angle(k2),polarization,thickness,nk_layer(k2,k1),nk_layer(k2+1,k1),k2,layers);
            if(k2==1)          
                tmm_matrix=coherent_matrix;
            else
                tmm_matrix=tmm_matrix*coherent_matrix;        
            end
        end
        r_coeff=tmm_matrix(2,1)/tmm_matrix(1,1);    t_coeff=1/tmm_matrix(1,1);
        if(polarization==0)
            Tran(k1)=abs(t_coeff)^2*real(conj(nk_layer(layers,1))*cosd(angle(layers)))/(nk_layer(1,1)*cosd(angle(1))); 
        elseif(polarization==90)
            Tran(k1)=abs(t_coeff)^2*real(nk_layer(layers,1)*cosd(angle(layers)))/(nk_layer(1,1)*cosd(angle(1)));   
        end
        Ref(k1)=abs(r_coeff)^2;
    end
end

function [coherent_matrix,angle1]=coherent_tmm(wv,angle0,polarization,thickness,n0,n1,layer_no,layers)
    [interface,angle1]=coherent_interface(angle0,polarization,n0,n1);
    if(layer_no<layers-1)
        [propagation]=coherent_propagation(wv,thickness(layer_no),n1,angle1);
        coherent_matrix=interface*propagation;
    else
        coherent_matrix=interface;
    end
end

function [interface,angle1]=coherent_interface(angle0,polarization,n0,n1)
    angle1=asind((n0/n1)*sind(angle0));
    [r,t]=Fresnel(polarization,[angle0 angle1],[n0 n1]);
    interface=(1/t)*[1 r; r 1];
end

function [propagation]=coherent_propagation(wv,thickness,n1,angle1)
    % incoherent/mixed coherent-incoherent layers formalism in random phase
    if(thickness<1)     
        delta=2*pi*(thickness/wv)*n1*cosd(angle1);
    else
        delta=2*pi*(thickness/wv)*n1*cosd(angle1)+(2*rand-1)*pi;
    end
    propagation=[exp(-1i*delta) 0; 0 exp(1i*delta)];
end

function [r,t]=Fresnel(polarization,angle,index)
    if(polarization==0)
        r=(index(1)*cosd(angle(2))-index(2)*cosd(angle(1)))./(index(1)*cosd(angle(2))+index(2)*cosd(angle(1)));
        t=2*index(1)*cosd(angle(1))./(index(1)*cosd(angle(2))+index(2)*cosd(angle(1)));
    elseif(polarization==90)
        r=(index(1)*cosd(angle(1))-index(2)*cosd(angle(2)))./(index(1)*cosd(angle(1))+index(2)*cosd(angle(2)));
        t=2*index(1)*cosd(angle(1))./(index(1)*cosd(angle(1))+index(2)*cosd(angle(2)));
    end
end