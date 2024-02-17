import React from "react";
import ReactApexChart from "react-apexcharts";
import { ChartDataPoint, URLParams } from "src/types";

interface SimpleTermChartProps {
  paramDataValue: {
    params: URLParams;
    data: ChartDataPoint[];
  };
}

const SimpleTermChart: React.FC<SimpleTermChartProps> = ({
  paramDataValue,
}) => {
  const chartData = {
    series: [
      {
        name: "Reports",
        data: paramDataValue.data,
      },
    ],
  };

  return (
    <div>
      <h3>{paramDataValue.params.sex}</h3>
      <h4>
        {paramDataValue.params.age?.min}-{paramDataValue.params.age?.max}
      </h4>
      <ReactApexChart
        options={{
          chart: {
            type: "bar",
          },
          xaxis: {
            categories: paramDataValue.data.map((d) => d.x),
          },
        }}
        type="bar"
        title="Term Data"
        series={chartData.series}
      />
    </div>
  );
};

export default SimpleTermChart;
