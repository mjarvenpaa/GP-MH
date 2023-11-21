function [] = write_err(id_seed, id_model, err)
% Saves the error that was observed during cluster run and the 
% corresponding settings to file. 

fid = fopen('errorfile.txt','a+');
fprintf(fid, '--------------------\n');
fprintf(fid, 'id_seed: %i, id_model: %i', id_seed, id_model);
fprintf(fid, '\n');
fprintf(fid, '%s', err.getReport('extended','hyperlinks','off'));
fprintf(fid, '\n--------------------\n');
fprintf(fid, '\n');
fclose(fid);
end
