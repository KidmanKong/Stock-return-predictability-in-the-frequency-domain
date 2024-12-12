function data_comp=wavelet_decomposing_function(data,wname,N)

% data��ԭʼ���ݣ�����ΪT������
% wname��С����
% N��С���任�ļ���

T=length(data);
data_comp=nan(T,N+1);
for i=1:N+1
    [C,L]=wavedec(data,N,wname);
    if i==1
        C(L(1)+1:end)=C(L(1)+1:end).*0;
    else
        C(1:sum(L(1:i))-sum(L(1:i-1)))=C(1:sum(L(1:i))-sum(L(1:i-1))).*0;
        C(sum(L(1:i))+1:end)=C(sum(L(1:i))+1:end).*0;
    end
    data_comp(:,i)=waverec(C,L,wname);
end