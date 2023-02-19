function [] = make_solnp()


platform = convertCharsToStrings(computer('arch'));

mexfname = 'SOLNP_plus';

% change your mkl path here
mkl_path = 'path to your mkl root';

% change your osqp path here
osqp_path = 'path to your osqp root';


if platform == "maci64"
    mkl_lib_path = fullfile(mkl_path, 'lib');
else
    mkl_lib_path = fullfile(mkl_path, 'lib', 'intel64');
end
osqp_lib_path = fullfile(osqp_path, 'build', 'out');



if exist(mkl_path) && exist(osqp_path)
    lib_path = [convertCharsToStrings(mkl_lib_path), convertCharsToStrings(osqp_lib_path)];
    lib_path = join('-L' + lib_path);
else
    error('Please add your pathes to mkl and osqp\n');
end 


link = '';

if platform == "win64"
    fprintf('Linking MKL in Windows \n');
    link = [link, ' -lmkl_intel_lp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "maci64"
    fprintf('Linking MKL in MacOS \n');
    link = [link, join(mkl_lib_path + ["/libmkl_intel_lp64.a","/libmkl_core.a", "/libmkl_sequential.a "])];
elseif platform == "glnxa64"
    fprintf('Linking MKL in Linux \n');

    sys_cmd = ['gcc -fPIC -I', matlabroot, '/extern/include ' , '-I./include -I', mkl_path, '/include -I', osqp_path, '/include ',...
                '-shared -o ', mexfname, '.mexa64 ', './source/*.c ', '-Wl,--start-group ', mkl_lib_path, '/libmkl_core.a ',...
                mkl_lib_path, '/libmkl_intel_lp64.a ', mkl_lib_path, '/libmkl_sequential.a ', mkl_lib_path, '/libmkl_intel_thread.a ',...
                osqp_lib_path, '/libosqp.a -Wl,--end-group -L', matlabroot, '/bin/glnxa64 -ldl -lpthread -lmx -lmex -lmat -lm'];

    system(sys_cmd);
    return
else
    error('Unsupported platform.\n');
end

link = [link, '-losqp'];

psrc = fullfile('.', 'source');  


src_files = dir( [psrc '/*.c'] );
srclist = [];

for i = 1:length(src_files)
    srclist = [srclist,convertCharsToStrings(src_files(i).name)];
end

src = fullfile(psrc, srclist);

pinc = "include";


mkl_include = fullfile(mkl_path, "include");
osqp_include = fullfile(osqp_path, "include");

inc = [pinc,mkl_include,osqp_include];
inc = join('-I' + inc);

cmd = 'mex -output ' + join([mexfname, lib_path, src, inc, link]);

fprintf('%s\n',cmd);
eval(replace(cmd, 'Program Files (x86)', "'Program Files (x86)'"));

end