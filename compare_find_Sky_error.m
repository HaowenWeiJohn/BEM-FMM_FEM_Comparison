
file_laptop = 'find_Sky_error_save3laptop.mat';
file_Sky = 'find_Sky_error_save3Sky.mat';

l = load(file_laptop);
s = load(file_Sky);

fn_l = fieldnames(l);
fn_s = fieldnames(s);
assert(numel(fn_l)==numel(fn_s));
for k=1:numel(fn_l)
    disp(fn_s{k});
    disp(fn_l{k});
    
    if (isnumeric(l.(fn_l{k})))
        temp = l.(fn_l{k}) - s.(fn_l{k});
        disp(norm(temp(:)));
    end
    
    if (isstruct(l.(fn_l{k})))
        fn_obj_l = fieldnames(l.(fn_l{k}));
        fn_obj_s = fieldnames(s.(fn_l{k}));
        for j=1:numel(fn_obj_l)
            disp(fn_obj_l{j});
            disp(fn_obj_s{j});
            if (isnumeric(l.(fn_l{k}).(fn_obj_l{j})))
                temp = l.(fn_l{k}).(fn_obj_l{j}) - s.(fn_l{k}).(fn_obj_l{j});
                disp(norm(temp(:)));
            end
        end
    end
    
    if (iscell(l.(fn_l{k})))
        for j=1:numel(l.(fn_l{k}))
            if (isnumeric(l.(fn_l{k}){j}))
                temp = l.(fn_l{k}){j} - l.(fn_l{k}){j};
                disp(norm(temp(:)));
            else
                disp(l.(fn_l{k}){j});
                disp(l.(fn_l{k}){j});
            end
        end
    end
end