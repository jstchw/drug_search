import { ApexOptions } from "apexcharts";
import { ThemeType } from "../../types";
import _ from "lodash";

const getCommonOptions = (theme: ThemeType) => {
  return {
    theme: {
      mode: theme as ThemeType,
    },
    chart: {
      toolbar: {
        show: false,
        tools: {
          zoom: false,
          zoomin: false,
          zoomout: false,
        },
      },
      background: theme === "dark" ? "#212529" : "",
    },
    yaxis: {
      labels: {
        style: {
          fontSize: "14px",
          colors: theme === "dark" ? "#ACB5BD" : "",
        },
      },
    },
  };
};

// const barColors = [
//   "#59768A",
//   "#035363",
//   "#32B2BF",
//   "#D5E0BE",
//   "#CE9062",
//   "#E0AB86",
//   "#C7CE8A",
//   "#6EB585",
//   "#325951",
//   "#6F9F9D",
// ];

export const getTermCarouselOptions = (theme: ThemeType): ApexOptions => {
  const commonOptions = getCommonOptions(theme);

  const specificOptions = {
    legend: {
      show: false,
    },
    chart: {
      type: "bar",
      width: "100%",
      toolbar: {
        show: false,
      },
      dropShadow: {
        enabled: true,
        top: 1,
        left: 1,
        blur: 2,
        color: theme === "dark" ? "#000" : "#000",
        opacity: 0.2,
      },
    },
    plotOptions: {
      bar: {
        horizontal: true,
        distributed: true,
        barHeight: "80%",
        borderRadius: 4,
        borderRadiusApplication: "end",
      },
    },
    grid: {
      show: false,
    },
    yaxis: {
      labels: {
        show: true,
        style: {
          fontSize: "14px",
          colors: theme === "dark" ? "#ACB5BD" : "",
        },
        formatter: (value: number) => {
          const stringValue = String(value);
          return (
            stringValue.charAt(0).toUpperCase() +
            stringValue.slice(1).toLowerCase()
          );
        },
      },
    },
    xaxis: {
      labels: {
        show: false,
      },
    },
  };

  return _.merge({}, commonOptions, specificOptions) as ApexOptions;
};

export const getTimeSeriesOptions = (theme: ThemeType): ApexOptions => {
  const commonOptions = getCommonOptions(theme);

  const specificOptions = {
    colors: ["#59768A"],
    legend: {
      show: true,
      position: "bottom",
    },
    chart: {
      type: "area",
    },
    stroke: {
      curve: "smooth",
    },
    fill: {
      type: "gradient",
    },
    dataLabels: {
      enabled: false,
    },
    xaxis: {
      type: "category",
      labels: {
        style: {
          fontSize: "14px",
          colors: theme === "dark" ? "#ACB5BD" : "",
        },
      },
      tooltip: {
        enabled: false,
      },
    },
    yaxis: {
      labels: {
        style: {
          fontSize: "14px",
          colors: theme === "dark" ? "#ACB5BD" : "",
        },
        formatter: (value: number) => {
          return value.toLocaleString();
        },
      },
    },
  };

  return _.merge({}, commonOptions, specificOptions) as ApexOptions;
};
