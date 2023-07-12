
import React from "react";
import {Button} from "react-bootstrap";
import './ReadMore.css';

export const ReadMore = ({children}) => {
    const [isExpanded, setIsExpanded] = React.useState(true);
    const isTooLong = children.length > 20;

    const buttonText = isExpanded ? "show more" : "show less";
    return (
        <React.Fragment>
            <p>
                {isExpanded ? children.slice(0, 20) : children}
            </p>
            <div className={'d-flex justify-content-center'}>
                {isTooLong &&
                    <Button
                        size={'sm'}
                        className={'expand-button'}
                        onClick={() => setIsExpanded(!isExpanded)}>{buttonText}
                    </Button>
                }
            </div>
        </React.Fragment>
    )
}