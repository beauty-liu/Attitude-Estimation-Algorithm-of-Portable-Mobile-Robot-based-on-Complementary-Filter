% classical_complementary_filters
classdef classical_complementary_filters < handle
%MAYHONYAHRS Madgwick's implementation of Mayhony's AHRS algorithm
%
%   For more information see:
%   http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
%
%   Date          Author          Notes
%   28/09/2011    SOH Madgwick    Initial release
 
    %% Public properties
    properties (Access = public)
        tao = 1;                     % algorithm proportional gain
        Ts = 0.01;                     % algorithm integral gain
        angle = 0; 
    end
    
    %% Public properties
    properties (Access = private)
        last_angle = 0;          
    end    
 
    %% Public methods
    methods (Access = public)
        function obj = classical_complementary_filters(varargin)
            for i = 1:2:nargin
                if  strcmp(varargin{i}, 'tao'), obj.tao = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Ts'), obj.Ts = varargin{i+1};
                else error('Invalid argument');
                end
            end;
        end

        function obj = UpdateIMU(obj, w_m,angle_m )
           obj.angle =  (obj.tao/(obj.tao + obj.Ts)) * ( obj.last_angle + obj.Ts*w_m) + (obj.Ts/(obj.tao + obj.Ts)) * angle_m;
           obj.last_angle = obj.angle;
        end
    end
end