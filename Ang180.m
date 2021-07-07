function [ang180] = Ang180(ang360)
if ang360>=180*pi/180
    while ang360>=180*pi/180
    ang360=ang360-360*pi/180;
    end
    ang180=ang360;
    return
end

if ang360<-180*pi/180
    while ang360<-180*pi/180
    ang360=ang360+360*pi/180;
    end
    ang180=ang360;
    return
end

ang180=ang360;
return