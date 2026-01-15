%% iman_statusreport(discordName, bigMood)
% Long story short, I had a bot that helped me try to get a GPU when it was hard...
% This is a chunk of that code, but used for good!
% It will ping the webhookUrl (a discord channel), tagging whoever the person that is 
% the "discordName" and change the message depending on "bigMood"
% 
% discordName - person's name / the channel / or the role as a string or cell array of names to ping
% example: 'nick' - that's me!
% current options....
%
% people: (in alphabetical order)
% christi - 889986108526575626;
% cua or cuauhtemoc -1090695735630102558
% daniel - 528284146867503124
% elijah - 172499874447753217
% florine - 1090738615887335554
% jason - 240689889258110976
% john - 1090430409491357807
% madhura - 1090425961889157131
% marion - 486640365819133953
% markhus - 1090705907849113734
% mike - 1091057639124901999
% nick - 257430809181421568
%
%
% roles (sub groups):
% metabolism - 1090520026592837692
% epi - 1090523156768358493 
% rnaseq - 1090534095911669812
% egfr - 1090706547526615171
%
% 
% bigMood - string with the mood of the message you want to send.
% example - 'happy celltracer'
% Current options for bigMood:
% (format: bigMood - what that message will say)
% sad celltracer - @discordName celltracer_v2 has crashed! :(
% happy celltracer - @disordName celltracer_v2 is all finished!
% meetign time - @discordName hey everyone its meeting time!
%
% Feel free to add more or steal the code, I guess?
%
% With love,
% Nick
% 2023
%
function iman_statusreport(discordName,bigMood)

% this is the webhook URL for our lab's discord channel aaybe fix this to make it so you can change this?
labDiscordWebhook = 'https://discord.com/api/webhooks/1093417323404202004/Ux8mdQlzYyRThP-iNCWWNB-60mtI9cYDcYxpwUFu-tFgYOpnpO1gRTMC30KkcnI3xMAn'; % discord webhook.... might be dangerous to have everywhere
options = weboptions('RequestMethod', 'post', 'MediaType', 'application/json'); % set options for sending a message via discord


% check if discordName is a string or a cell array of strings
if ischar(discordName) || isstring(discordName)
   discordName = {discordName};
elseif iscellstr(discordName)
else; return
end

discordList = ''; % here is the list of names

for iName = 1:numel(discordName)
    switch lower(discordName{iName})
        % people (alphabetical order)
        case 'christi'; nameToAdd = '889986108526575626';
        case {'cua', 'cuauhtemoc'}; nameToAdd = '1090695735630102558';
        case 'daniel';    nameToAdd = '528284146867503124';
        case 'elijah';    nameToAdd = '172499874447753217';
        case 'florine';   nameToAdd = '1090738615887335554';
        case 'jason';     nameToAdd = '240689889258110976';
        case 'john';      nameToAdd = '1090430409491357807';
        case 'madhura';   nameToAdd = '1090425961889157131';
        case 'marion';    nameToAdd = '486640365819133953';     
        case 'markhus';   nameToAdd = '1090705907849113734';
        case 'mike';      nameToAdd = '1091057639124901999';
        case 'nick';      nameToAdd = '257430809181421568';

        % subgroups
        case 'metabolism'; nameToAdd = '1090520026592837692';
        case 'epi';        nameToAdd = '1090523156768358493';
        case 'rnaseq';     nameToAdd = '1090534095911669812';
        case 'egfr';       nameToAdd = '1090706547526615171';

        otherwise % get out of there
            continue
    end
    discordList = [discordList,'<@', nameToAdd, '> ']; % add them to the ping
    
end

if isempty(discordList); return; end % don't send an empty ping

switch bigMood
    case 'sad celltracer' % if your celltracer_v2 broke
        message = [discordList, 'celltracer_v2 has crashed! :('];
    case 'happy celltracer' % if your celltracer_v2 did it!
        message = [discordList, 'celltracer_v2 is all finished!'];
    case 'meeting time'
        message = [discordList, 'hey everyone its time to meet!'];
    otherwise % idk why youre calling me, but here we are
        message = [discordList, 'idk friend, you brought me here and now we are in a weird position because you didnt tell me why...'];
end
message = struct('content',message); % convert it to a struct
webwrite(labDiscordWebhook, message, options); % SEND IT

end % iman_statusreport