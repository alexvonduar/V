function centered_coordinates = coordinates_center(coordinates, width, height, focal)
center = [width, height] / 2;
if size(coordinates, 2) == 2
  centered_coordinates = bsxfun(@minus, coordinates, center) / focal;
else
  centered_coordinates = bsxfun(@minus, coordinates, [center, center]) / focal;
end