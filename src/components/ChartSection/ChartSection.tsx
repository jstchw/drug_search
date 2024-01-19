import CasesTimeSeriesChart from "../ApexChart/CasesTimeSeriesChart"
import {OverlayTrigger, Popover} from "react-bootstrap"
// import WordCloudChart from "../WordCloudChart/WordCloudChart"
import { ExclamationCircle, QuestionCircle } from "react-bootstrap-icons"
import TermChart from "../ApexChart/TermChart";

const getChartWarning = (linkText: string, text: string, status: string) => {
    const popover = (
        <Popover>
            <Popover.Body>
                <div>{text}</div>
            </Popover.Body>
        </Popover>
    )

    return (
        <OverlayTrigger trigger={['hover', 'focus']} placement={'bottom'} overlay={popover}>
            <div style={{cursor: 'default'}}>
                <span>{linkText}</span>
                <sup>
                    {status === 'warn' &&
                        <ExclamationCircle className={'mx-1'}/>}
                    {status === 'info' &&
                        <QuestionCircle className={'mx-1'}/>}

                </sup>
            </div>
        </OverlayTrigger>
    )
}

const ChartSection = () => {
    const adverse_effects_text = 'Correlation does not imply causation.'

    return (
        <div className={'mt-4'}>
            <div>
                <h4 className={'text-center'}>Reports over time</h4>
                <CasesTimeSeriesChart />
            </div>
            <div className={'mt-4'}>
                <h4 className={'text-center'}>
                    {getChartWarning('Adverse effects', adverse_effects_text, 'warn')}
                </h4>
                <TermChart />
            </div>
            {/*<div className={'mt-4'}>*/}
            {/*    <h4>Word Cloud</h4>*/}
            {/*    <WordCloudChart words={props.searchResults} searchOptions={props.searchOptions}/>*/}
            {/*</div>*/}
        </div>
    );
}

export default ChartSection;