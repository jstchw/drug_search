import React from 'react';
import ReactWordcloud from 'react-wordcloud';
import 'tippy.js/dist/tippy.css';
import 'tippy.js/animations/scale.css';

const options = {
    enableTooltip: true,
    enableOptimizations: true,
    deterministic: true,
    fontFamily: "trebuchet ms",
    fontSizes: [5, 60],
    fontStyle: "normal",
    fontWeight: "normal",
    padding: 1,
    rotations: 1,
    rotationAngles: [0, 90],
    scale: "sqrt",
    spiral: "rectangular",
    transitionDuration: 1000
};

const WordCloudChart = (props) => {
    const data = Object.entries(props.ADEArray)
        .map(([text, value]) => ({text, value}))
        .slice(0, 50)

    const callbacks = {
        getWordColor: word => {
            const maxFrequency = Math.max(...data.map(w => w.value)); // get the maximum frequency
            const frequencyRatio = word.value / maxFrequency; // calculate the ratio of the word's frequency to the maximum frequency
            const redComponent = Math.floor(frequencyRatio * 255); // calculate the red component (0-255)
            return `rgb(${redComponent}, 0, 0)`; // return an RGB color string
        },
        // You can add more callbacks if you want...
    };

    return (
        <ReactWordcloud words={data} options={options} callbacks={callbacks} />
    )
}

export default WordCloudChart;