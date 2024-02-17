import React from "react";
import ReactApexChart from "react-apexcharts";
import { ChartDataPoint, ThemeType, URLParams } from "src/types";
import { getTermCarouselOptions } from "./chartOptions";
import { ThemeContext } from "../../contexts/ThemeContext";

const formatParamsForDisplay = (params: URLParams) => {
  return {
    ...params,
    sex: params.sex && params.sex.charAt(0).toUpperCase() + params.sex.slice(1),
    age: params.age && `${params.age.min}-${params.age.max}`,
  };
};

interface SimpleTermChartProps {
  paramDataValue: {
    params: URLParams;
    data: ChartDataPoint[];
  };
}

const SimpleTermChart: React.FC<SimpleTermChartProps> = ({
  paramDataValue,
}) => {
  const { theme } = React.useContext(ThemeContext);

  const formattedParams = formatParamsForDisplay(paramDataValue.params);

  const chartData = {
    series: [
      {
        name: "Reports",
        data: paramDataValue.data,
      },
    ],
  };

  const simpleChartOptions = getTermCarouselOptions(theme as ThemeType);

  return (
    <div className={"text-center"}>
      <h3>{formattedParams.sex}</h3>
      <h5>{formattedParams.age}</h5>
      <ReactApexChart
        options={simpleChartOptions}
        type="bar"
        title="Term Data"
        series={chartData.series}
      />
    </div>
  );
};

export default SimpleTermChart;
