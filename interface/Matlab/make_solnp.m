function [] = make_solnp(LINK_MKL)

mexfname = 'SOLNP_plus';

% change your mkl path here
mkl_path = "YOUR MKL PATH";
% change your lapack path here
lapack_path = "Your LAPACK path";
% change your osqp path here
osqp_path = "YOUR OSQP PATH";

mkl_lib_path = fullfile(mkl_path, 'lib', 'intel64');
osqp_lib_path = fullfile(osqp_path, 'lib');
lapack_lib_path = fullfile(lapack_path, 'lib');

if exist(mkl_path) && exist(osqp_path)
    if LINK_MKL
        lib_path = [convertCharsToStrings(mkl_lib_path), convertCharsToStrings(osqp_lib_path)];
    else
        lib_path = [convertCharsToStrings(lapack_lib_path), convertCharsToStrings(osqp_lib_path)];
    end
    lib_path = join('-L' + lib_path);
    platform = convertCharsToStrings(computer('arch'));
else
    error('Please add your pathes to osqp and mkl or lapack\n');
end 


link = '';

if platform == "win64"
    if LINK_MKL
        fprintf('Linking MKL on Windows \n');
        link = [link, ' -lmkl_intel_lp64',...
            ' -lmkl_core', ' -lmkl_sequential '];
    else
        fprintf('Linking LAPACK on Windows \n');
        link = [link, ' -lblas', ' -llapack '];
    end
elseif platform == "maci64"
    if LINK_MKL
        fprintf('Linking MKL on MacOS \n');
        link = [link, ' -lmkl_intel_lp64',...
            ' -lmkl_core', ' -lmkl_sequential '];
    else
        fprintf('Linking LAPACK on MacOS \n');
        link = [link, ' -lblas', ' -llapack '];
    end
    
elseif platform == "glnxa64"
    if LINK_MKL
        fprintf('Linking MKL in Linux \n');
        sys_cmd = ['gcc -fPIC -I', matlabroot, '/extern/include ' , '-I../../include -I', mkl_path, '/include -I', osqp_path, '/include ',...
                    '-shared -o ', mexfname, '.mexa64 ', '../../source/*.c ', '-Wl,--start-group ', mkl_lib_path, '/libmkl_core.a ',...
                    mkl_lib_path, '/libmkl_intel_lp64.a ', mkl_lib_path, '/libmkl_sequential.a ', mkl_lib_path, '/libmkl_intel_thread.a ',...
                    osqp_lib_path, '/libosqp.a -Wl,--end-group -L', matlabroot, '/bin/glnxa64 -ldl -lpthread -lmx -lmex -lmat -lm'];
    else
        fprintf('Linking LAPACK in Linux \n');
        sys_cmd = ['gcc -fPIC -I', matlabroot, '/extern/include ' , '-I../../include -I', osqp_path, '/include ','-shared -o ', mexfname, ...
                    '.mexa64 ', '../../source/*.c ', '-Wl,--start-group ', lapack_lib_path, '/libblas.a ', lapack_lib_path, '/liblapack.a ', ...
                    osqp_lib_path, '/libosqp.a -Wl,--end-group -L', matlabroot, '/bin/glnxa64 -ldl -lpthread -lmx -lmex -lmat -lm'];      
    end
    system(sys_cmd);
    return
else
    error('Unsupported platform.\n');
end

link = [link, '-losqp'];

psrc = fullfile('..', '..', 'source');  

src_files = dir( [psrc '\*.c'] );
srclist = [];

for i = 1:length(src_files)
    srclist = [srclist,convertCharsToStrings(src_files(i).name)];
end

src = fullfile(psrc, srclist);

pinc = fullfile('..', '..', 'include');
osqp_include = osqp_path + "\include";

inc = [pinc, osqp_include];
inc = join('-I' + inc);

cmd = 'mex -g -output ' + join([mexfname, lib_path, src, inc, link]);

fprintf('%s\n',cmd);
eval(replace(cmd, 'Program Files (x86)', "'Program Files (x86)'"));

end