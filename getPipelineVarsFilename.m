% Small helper function which returns the enclosing directory of the
% pipeline code by splitting the filename off mfilename('fullpath')
function fname = getPipelineVarsFilename()
    varpath = mfilename('fullpath');
    splitpath = split(varpath,'\');
    name = join([ join(splitpath(1:length(splitpath)-1), '\') 'shared_vars.mat'], '\');
    fname = name{1};
end