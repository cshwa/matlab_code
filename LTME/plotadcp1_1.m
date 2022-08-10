% plot the adcp data as direction arrows.  Loads the .mat files created by
% adcpdata2array.m.  Must be run from the same directory the files are in.
% Bill Scuba   Oct. 15, 2003 and edit by Peter(2004-4-23)
% INPUT data: 현재 디렉토리에 있는 모든 *.mat 파일
% m-file: plotadcp1.m
% array별로 저장된 *.mat파일을 읽어 배열 ping의 데이터를 해류벡터로 도시
% 먼저 adcpdata2array.m을 실행 후 본 실행파일을 실행하면 수평해류분포도가 작성
% 

clear dirstruct            % 배열 dirstruct에 저장된 내용을 삭제

figure                     % 그림 창 (1개)을 연다
hold on                    % 모든 그림을 동일한 창에 그린다
dirstruct=dir('*.mat')     % 현재 디렉토리의 확장자가 mat인 모든 파일을 배열 dirstruct에 저장

% 현재 디렉토리의 확장자가 mat인 모든 파일을 반복 처리하는 루틴
	for yi=1:length(dirstruct)   % yi는 1에서 배열 dirstruct의 파일갯수 번호까지
    clear ping                       %배열 ping에 기억된 내용을 지움 
    load(dirstruct(yi).name) 
	for xi=1:50:length(ping)     % xi는 1행부터 ping의 전체행 번호까지 50번째 행 번호마다 처리하는 루프
		quiver(ping(xi).lon,ping(xi).lat,ping(xi).adcpdata(1,4),ping(xi).adcpdata(1,5),.002)
		% 배열 ping의 xi번째 방의 이름인 배열 lon (x값: 경도값), 
		% 배열 ping의 xi번째 방의 이름인 배열 lat (y값: 위도값),
		% 배열 ping의 xi번째 방의 이름인 배열 adcpdata의 1행, 4열의 값 (동방성분)
		% 배열 ping의 xi번째 방의 이름인 배열 adcpdata의 1행, 5열의 값 (북방성분))
		% 화살표의 길이 (length of arrow)
	end
    fprintf('u = %g \n',sqrt(ping(xi).adcpdata(1,4)^2 +ping(xi).adcpdata(1,5)^2))
		% u 유속 값은 동방성분(U)과 북방성분(V)의 각각의 제곱을 합한값에 제곱근을 통해 구함
    end
hold off           % 모든 그림을 동일한 창에 그리지 않음