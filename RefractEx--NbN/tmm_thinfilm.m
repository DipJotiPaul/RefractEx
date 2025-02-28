function [spectra]=tmm_thinfilm(wv,incidence_angle,polarization,thickness,nk_layer,nk_substrate,spectra_no)
    if(max(thickness)<1)                                 % coherent layers formalism when thickness is below 1um
        if(polarization==0 || polarization==90) 
            spectra=tmm(wv,incidence_angle,polarization,thickness,nk_layer,nk_substrate,spectra_no);
        else
            spectra1=tmm(wv,incidence_angle,0,thickness,nk_layer,nk_substrate,spectra_no);
            spectra2=tmm(wv,incidence_angle,90,thickness,nk_layer,nk_substrate,spectra_no);
            spectra=0.5*(spectra1+spectra2);
        end
    else                                                               % incoherent layers formalism in random phase when thickness > 1um
        iteration=50;           temp=0;
        for k1=1:iteration
            if(polarization==0 || polarization==90)
                 spectra=tmm(wv,incidence_angle,polarization,thickness,nk_layer,nk_substrate,spectra_no);
            else
                spectra1=tmm(wv,incidence_angle,0,thickness,nk_layer,nk_substrate,spectra_no);
                spectra2=tmm(wv,incidence_angle,90,thickness,nk_layer,nk_substrate,spectra_no);
                spectra=0.5*(spectra1+spectra2);
            end
        temp=temp+spectra;    
        end
        spectra=temp/iteration;
    end
    spectra=spectra*100;
end

function [spectra]=tmm(wv,incidence_angle,polarization,thickness,nk,nk_substrate,spectra_no)
    index1=1;           index2=nk;          index3=nk_substrate;    index4=1;
    angle1=incidence_angle;                            angle2=asind((index1/index2)*sind(angle1));
    if(polarization==0)
        t=2*index1*cosd(angle1)./(index1*cosd(angle2)+index2*cosd(angle1));
        r=(index1*cosd(angle2)-index2*cosd(angle1))./(index1*cosd(angle2)+index2*cosd(angle1));
    elseif(polarization==90)
        t=2*index1*cosd(angle1)./(index1*cosd(angle1)+index2*cosd(angle2));
        r=(index1*cosd(angle1)-index2*cosd(angle2))./(index1*cosd(angle1)+index2*cosd(angle2));
    end
    interface1=(1/t)*[1 r; r 1];
    if(thickness(1)<1)
        delta=2*pi*(thickness(1)/wv)*index2*cosd(angle2);
    else
        delta=2*pi*(thickness(1)/wv)*index2*cosd(angle2)+(2*rand-1)*pi;
    end
    propagation1=[exp(-1i*delta) 0; 0 exp(1i*delta)];
    angle3=asind((index2/index3)*sind(angle2));
    if(polarization==0)
        t=2*index2*cosd(angle2)./(index2*cosd(angle3)+index3*cosd(angle2));
        r=(index2*cosd(angle3)-index3*cosd(angle2))./(index2*cosd(angle3)+index3*cosd(angle2));
    elseif(polarization==90)
        t=2*index2*cosd(angle2)./(index2*cosd(angle2)+index3*cosd(angle3));
        r=(index2*cosd(angle2)-index3*cosd(angle3))./(index2*cosd(angle2)+index3*cosd(angle3));
    end
    interface2=(1/t)*[1 r; r 1];
    if(thickness(2)<1)
        delta=2*pi*(thickness(2)/wv)*index3*cosd(angle3);
    else
        delta=2*pi*(thickness(2)/wv)*index3*cosd(angle3)+(2*rand-1)*pi;
    end
    propagation2=[exp(-1i*delta) 0; 0 exp(1i*delta)];
    angle4=asind((index3/index4)*sind(angle3));
    if(polarization==0)
        t=2*index3*cosd(angle3)./(index3*cosd(angle4)+index4*cosd(angle3));
        r=(index3*cosd(angle4)-index4*cosd(angle3))./(index3*cosd(angle4)+index4*cosd(angle3));
    elseif(polarization==90)
        t=2*index3*cosd(angle3)./(index3*cosd(angle3)+index4*cosd(angle4));
        r=(index3*cosd(angle3)-index4*cosd(angle4))./(index3*cosd(angle3)+index4*cosd(angle4));
    end
    interface3=(1/t)*[1 r; r 1];
    tmm_matrix=interface1*propagation1*interface2*propagation2*interface3;
    r_coeff=tmm_matrix(2,1)/tmm_matrix(1,1);    t_coeff=1/tmm_matrix(1,1);
    if(polarization==0 && spectra_no==1)
        spectra=abs(t_coeff)^2*real(conj(index4)*cosd(angle4))/(index1*cosd(angle1)); 
    elseif(polarization==90 && spectra_no==1)
        spectra=abs(t_coeff)^2*real(index4*cosd(angle4))/(index1*cosd(angle1));   
    elseif(spectra_no==2)
        spectra=abs(r_coeff)^2;
    end
end