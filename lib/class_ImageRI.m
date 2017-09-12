classdef class_ImageRI
    properties
        im
        F
        change_pos
        O
        param
    end
    methods
    function im_ = getIm(obj, k)
        if obj.param.use_image_variation
            for o = 1:obj.O
                im_ = obj.im;
                im_patch = im_(obj.param.var_region_centre{o}(1) + obj.change_pos, ...
                    obj.param.var_region_centre{o}(2) + obj.change_pos);
                im_patch = im_patch .* obj.F;
                im_patch = obj.param.sin_var.max_amplitude * im_patch * ...
                    sin(obj.param.sin_var.phase{o} + k * 2 * pi / obj.param.sin_var.period{o});
                im_(obj.param.var_region_centre{o}(1) + obj.change_pos, ...
                    obj.param.var_region_centre{o}(2) + obj.change_pos) = im_(obj.param.var_region_centre{o}(1) + obj.change_pos, ...
                    obj.param.var_region_centre{o}(2) + obj.change_pos) + im_patch;
                im_(im_ < 0) = 0;
            end
        end
        if ~obj.param.use_image_variation
            im_ = obj.im;
        end
    end
    
    function obj = setIm(obj, im_)
        obj.im = im_;
    end
        
    function obj = setParams(obj, param)
        if strcmp(param.var_region_shape, 'gaussian')
            obj.param = param;
            obj.O = length(param.var_region_centre);
            sxy = param.var_region_size / (2*sqrt(2*log(2)));
            obj.change_pos = -param.var_region_size:1:param.var_region_size;
            [X1, X2] = meshgrid(obj.change_pos, obj.change_pos);

            theta = 0;

            a = cos(theta)^2/(2*sxy^2) + sin(theta)^2/(2*sxy^2);
            b = -sin(2*theta)/(4*sxy^2) + sin(2*theta)/(4*sxy^2);
            c = sin(theta)^2/(2*sxy^2) + cos(theta)^2/(2*sxy^2);

            obj.F = exp( - (a*X1.^2 - 2*b*X1.*X2 + c*X2.^2));
            obj.F = obj.F./norm(obj.F(:));
        end
    end
    end
end

