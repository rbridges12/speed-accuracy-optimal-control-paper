function time = read_time(file)

%file = 'ue_full_bestSoFar_states.sto';
f = fopen(file);
time = [];

for i = 1 : 6
    line = fgets(f);
end

while ~feof(f)
     line = fgets(f);
     line= strsplit(line);
     num = str2double(line{1});
     time = [time; num];
end