;; this function converts the decimals in a number to 'p', for use in file names
;; example: the input '1.2' is output as '1p2'

FUNCTION num_formatter_filename, input

  output = number_formatter(input)
  for i = 0, n_elements(input)-1 do begin
    number = output[i]
    position = strpos(number, '.')
    if position ne -1 then number = strmid(number,0,position) + 'p' + strmid(number, position+1, strlen(number)-position)
    output[i] = number
  endfor
  return, output
  
END