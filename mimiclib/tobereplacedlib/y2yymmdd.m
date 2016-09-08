function [YYMMDD]=y2yymmdd(date)

[yr,mm,dd] = jd2cal(yr2jd(date));

YYMMDD = [sprintf('%02.f',yr),sprintf('%02.0f',mm),sprintf('%02.0f',dd)];
YYMMDD = YYMMDD(3:end);

