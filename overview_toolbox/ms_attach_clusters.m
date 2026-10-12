function S = ms_attach_clusters(S, clustLatent, clustCyCIF, varargin)
%CT_ATTACH_CLUSTERS  Attach latent-trace and CyCIF cluster labels to S.cellT.
%
%   Adds two label columns so you can cross-tab "clusters by live-trace
%   dynamics" against "clusters by endpoint protein state" on the SAME cells.
%   Both are row-aligned to S (and thus to Lat/data/cycifT).
%
% USAGE:
%   % latent labels already row-aligned to S.Lat; cycif labels keyed by (xy,cellid):
%   S = ms_attach_clusters(S, clustID_latent, cycifClustTable);
%   % both already row-aligned vectors:
%   S = ms_attach_clusters(S, clustID_latent, clustID_cycif_aligned, 'CyCIFAligned',true);
%
% INPUT:
%   S            struct from dl_build_live_cycif_tsc
%   clustLatent  [nCells x 1] latent cluster labels, ROW-ALIGNED to S.Lat
%                (from ms_ae_cluster on S.Lat). Pass [] to skip.
%   clustCyCIF   EITHER a table with columns {rep(optional),xy,cellid,<label>}
%                to LEFT-JOIN, OR a [nCells x 1] row-aligned vector if
%                'CyCIFAligned' is true. Pass [] to skip.
%
% NAME-VALUE:
%   'CyCIFAligned'  clustCyCIF is already a row-aligned vector [false]
%   'CyCIFLabelCol' name of the label column in the cycif table ['cluster']
%   'CyCIFHasRep'   cycif table has a 'rep' column to match on [false]
%                   (false -> match on (xy,cellid) only; requires cycif
%                    labels came from a SINGLE rep or xy is globally unique)
%   'Verbose'       [true]
%
% OUTPUT:
%   S with S.cellT.clust_latent and/or S.cellT.clust_cycif added, plus
%   S.audit.clust_latent_n / .clust_cycif_matchrate for provenance.

q.CyCIFAligned=false; q.CyCIFLabelCol='cluster'; q.CyCIFHasRep=false; q.Verbose=true;
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'cycifaligned',  q.CyCIFAligned=logical(varargin{k+1}); k=k+2;
        case 'cyciflabelcol', q.CyCIFLabelCol=varargin{k+1};         k=k+2;
        case 'cycifhasrep',   q.CyCIFHasRep=logical(varargin{k+1});  k=k+2;
        case 'verbose',       q.Verbose=logical(varargin{k+1});      k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

n = height(S.cellT);

% ---- latent labels (must be row-aligned) -------------------------------
if ~isempty(clustLatent)
    assert(numel(clustLatent)==n, ...
        'clustLatent has %d rows but S has %d — must be row-aligned to S.Lat.', ...
        numel(clustLatent), n);
    S.cellT.clust_latent = clustLatent(:);
    S.audit.clust_latent_n = numel(unique(clustLatent));
    if q.Verbose
        fprintf('attached clust_latent: %d clusters over %d cells.\n', ...
            S.audit.clust_latent_n, n);
    end
end

% ---- CyCIF labels (aligned vector OR keyed table to left-join) ---------
if ~isempty(clustCyCIF)
    if q.CyCIFAligned
        assert(numel(clustCyCIF)==n, 'aligned clustCyCIF must be [nCells x 1].');
        S.cellT.clust_cycif = clustCyCIF(:);
        S.audit.clust_cycif_matchrate = 1;
    else
        T = clustCyCIF;
        hasUid = ismember('uid', T.Properties.VariableNames);
        if hasUid
            req = {'uid', q.CyCIFLabelCol};
        else
            req = {'xy','cellid', q.CyCIFLabelCol};
        end
        assert(all(ismember(req, T.Properties.VariableNames)), ...
            'cycif cluster table needs columns: %s', strjoin(req,', '));

        % Preferred join: uid string (as authored by cycif_pool_replicates,
        % e.g. "rep1_2_1" = rep_xy_cellid). S has numeric rep/origXY/origCellID,
        % so reconstruct the SAME uid and match strings — avoids fragile numeric
        % coercion of prefixed reps ("rep1") and mixed types.
        if ismember('uid', T.Properties.VariableNames)
            keyT = string(T.uid);
            keyS = "rep" + string(S.cellT.rep) + "_" + ...
                   string(S.cellT.origXY) + "_" + string(S.cellT.origCellID);
        elseif q.CyCIFHasRep
            assert(ismember('rep',T.Properties.VariableNames), ...
                'CyCIFHasRep=true but table has no rep column.');
            keyT = local_num(T.rep)*1e9 + local_num(T.xy)*1e7 + local_num(T.cellid);
            keyS = local_num(S.cellT.rep)*1e9 + local_num(S.cellT.origXY)*1e7 + local_num(S.cellT.origCellID);
        else
            keyT = local_num(T.xy)*1e7 + local_num(T.cellid);
            keyS = local_num(S.cellT.origXY)*1e7 + local_num(S.cellT.origCellID);
            if numel(unique(S.cellT.rep))>1
                warning(['CyCIFHasRep=false but S has multiple reps — (xy,cellid) ' ...
                    'is NOT unique across reps. Labels may be mis-joined. Add a rep ' ...
                    'column and set CyCIFHasRep=true, or provide uid.']);
            end
        end

        [tf_, loc] = ismember(keyS, keyT);
        % cluster labels may be string ("1","noise") or numeric; keep as the
        % source type so "noise" survives. Store into a string column.
        srcLab = T.(q.CyCIFLabelCol);
        lab = strings(n,1); lab(:) = "";                 % unmatched -> "" (missing)
        lab(tf_) = string(srcLab(loc(tf_)));
        S.cellT.clust_cycif = lab;
        S.audit.clust_cycif_matchrate = mean(tf_);
        if q.Verbose
            fprintf('attached clust_cycif via join: %d/%d matched (%.1f%%), %d unmatched -> NaN.\n', ...
                sum(tf_), n, 100*mean(tf_), sum(~tf_));
        end
        if mean(tf_) < 0.9
            warning('clust_cycif match rate %.1f%% — many cells lack a CyCIF cluster label.', ...
                100*mean(tf_));
        end
    end
end

% ---- concordance preview (if both attached) ----------------------------
if all(ismember({'clust_latent','clust_cycif'}, S.cellT.Properties.VariableNames)) && q.Verbose
    cc = S.cellT.clust_cycif;
    if isstring(cc) || iscellstr(cc)
        ok = strlength(string(cc)) > 0;          % "" = unmatched/missing
    else
        ok = ~isnan(cc);
    end
    [ct,~,~,~] = crosstab(S.cellT.clust_latent(ok), cellstr(string(cc(ok))));
    fprintf('\nlatent x CyCIF cluster contingency (%d cells with both labels):\n', sum(ok));
    disp(ct);
    fprintf(['Read: strong diagonal/block structure = trace dynamics predict\n' ...
             'endpoint protein state. Uniform = the two labelings are independent.\n']);
end
end

% ---- helper: coerce numeric | numeric-string | "rep1"-style to double ----
function v = local_num(x)
    if isnumeric(x), v = double(x); return; end
    s = string(x);
    s = regexprep(s, '[^0-9\.\-]', '');   % strip non-numeric (e.g. "rep" in "rep1")
    v = str2double(s);
end
