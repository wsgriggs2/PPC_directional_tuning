function path = get_python_path

%% get user & computer names
user = getenv('username');
computer = getenv('computername');

% some older Unix OS versions will return empty arrays for getenv...
if isempty(user) || isempty(computer)
    try
        user = char(java.lang.System.getProperty('user.name'));
        computer = char(java.net.InetAddress.getLocalHost.getHostName);
        [~, mac_address] = system('ifconfig en0 | grep ether');
    catch
        error('Could not get username')
    end
end

%% set path according to user & computer

% Whitney     
if strcmp(user,'wsgriggs') && contains(computer, 'RHDVFAL')
        path = 'C:\Users\wsgriggs\AppData\Local\Programs\Python\Python311\python.exe';

% Whitney on Mac
elseif strcmp(user, 'wsgriggs') && contains(computer, 'Whitneys-MacBook')
   path = '/Users/wsgriggs/miniforge3/bin/python';
elseif strcmp(user, 'wsgriggs') && contains(mac_address, 'f8:4d:89:6f:ea:f6')
   path = '/Users/wsgriggs/miniforge3/bin/python';
end

% make sure our file separators are correct
path = strrep(path,'/',filesep);
path = strrep(path,'\',filesep);


end